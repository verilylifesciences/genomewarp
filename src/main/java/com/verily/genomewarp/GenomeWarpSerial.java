/*
 * Copyright 2016 Verily Life Sciences LLC.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.verily.genomewarp;

import static java.nio.charset.StandardCharsets.UTF_8;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.google.genomics.v1.Range;
import com.google.genomics.v1.Variant;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.RegionType;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.TargetStrand;
import com.verily.genomewarp.utils.Fasta;
import com.verily.genomewarp.utils.GenomeRange;
import com.verily.genomewarp.utils.GenomeRangeUtils;
import com.verily.genomewarp.utils.GenomeWarpUtils;
import com.verily.genomewarp.utils.GvcfToVcfAndBed;
import com.verily.genomewarp.utils.VariantToVcf;
import com.verily.genomewarp.utils.VcfToVariant;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderVersion;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.SortedMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Serial version of the genome warp workflow.
 *
 * Note: Due to a bug, HTSJDK cannot distinguish between PASS
 * and unfiltered for variant calls. As such, any filter fields
 * which are PASS are removed. Only filters which are non-PASS
 * remain in the output VCF.
 */
public final class GenomeWarpSerial {

  private static class Arguments {
    @Parameter(
      description = "Path to the chain file for liftover",
      names = "--lift_over_chain_path"
    )
    public String liftOverChainPath = null;

    @Parameter(
      description = "Path to uncompressed raw query VCF file",
      names = "--raw_query_vcf"
    )
    public String rawQueryVcf = null;

    @Parameter(
      description = "Path to uncompressed raw query BED file",
      names = "--raw_query_bed"
    )
    public String rawQueryBed = null;

    @Parameter(
      description = "Path to reference query fasta for BED file",
      names = "--ref_query_fasta"
    )
    public String refQueryFASTA = null;

    @Parameter(
      description = "Path to reference target fasta for BED file",
      names = "--ref_target_fasta"
    )
    public String refTargetFASTA = null;

    @Parameter(
      description = "Path to the working directory for intermediate and final outputs",
      names = "--work_dir"
    )
    public String workDir = null;

    @Parameter(
      description = "Keep homozygous reference calls",
      names = "--keep_homozygous_reference_calls"
    )
    public Boolean keepHomozygousReferenceCalls = false;

    @Parameter(
      description = "Path to output vcf holding the transformed variants",
      names = "--output_variants_file"
    )
    public String outputVariantsFile = null;

    @Parameter(
      description = "Path to output BED file holding the valid output regions",
      names = "--output_regions_file"
    )
    public String outputRegionsFile = null;

    @Parameter(description = "Remove the info field", names = "--remove_info_field")
    public Boolean removeInfoField = false;

    @Parameter(description = "Only run GenomeWarp. This options requires having the file "
      + "out.bed.annotated in the directory specified by the work_dir flag.",
      names = "--only_genome_warp")
    public Boolean onlyGenomeWarp = false;

    @Parameter(description = "VCF part chunk sizes", names = "--vcf_chunk_size")
    public Integer vcfChunkSize = Integer.valueOf(50000);

    @Parameter(description = "Generate intermediate files", names = "--intermediate_files")
    public Boolean intermediateFiles = false;

    @Parameter(description = "Min match parameter for lift over", names = "--min_match")
    public Float minMatch = 1.0f;

    @Parameter(description = "Target assembly short name", names = "--target_assembly")
    public String targetAssembly = "B38";

    @Parameter(description = "Species this VCF file belongs to", names = "--species")
    public String species = "Homo sapiens";

    @Parameter(description = "Path to uncompressed raw query gVCF file", names = "--raw_query_gvcf")
    public String rawQueryGvcf = null;

    @Parameter(description = "Window size for splitting BED file regions. "
        + "Smaller window size helps to avoid having regions with complex variations. It may improve the throughput",
        names = "--bed_window_size")
    public Integer bedWindowSize = 10000;

    @Parameter(description = "Skip improved regions preprocessing. If enabled, the output will be equivalent to v1.0.0. "
      + "Otherwise GenomeWarp splits bed regions in order to improve handling of complex genome transformations",
      names = "--simplified_regions_preproc")
    public Boolean simplifiedRegionsPreprocessing = false;
}

  // Used exclusively to facilitate storing
  private class RefNameToLength {
    public final String chr;

    public final Long length;

    public RefNameToLength(String chr, Long length) {
      this.chr = chr;
      this.length = length;
    }
  }

  public static final Logger logger = Logger.getLogger(GenomeWarpSerial.class.getName());

  private static final Arguments ARGS = new Arguments();

  private static final VCFHeaderVersion DEFAULT_VCF_VERSION = VCFHeaderVersion.VCF4_2;

  private static final String GENOME_WARP_VERSION = "GenomeWarp_v1.1.0";

  private static final String VCF_PART_NAME = "/vcfChunk.vcf.part";

  private static List<String> retrieveVcfHeader(String fileName) throws IOException {
    List<String> headerStrings = new ArrayList<>();
    String line;
    try (BufferedReader br = Files.newBufferedReader(Paths.get(fileName), UTF_8)) {
      while ((line = br.readLine()) != null) {
        if (line.length() > 0 && line.charAt(0) == '#') {
          headerStrings.add(line);
        } else {
          break;
        }
      }
    }

    return headerStrings;
  }


  public static RegionType getRegionType(HomologousRange range, Fasta queryFasta,
      Fasta targetFasta) {
    Range queryRange = range.getQueryRange();
    Range targetRange = range.getTargetRange();
    boolean isPositiveStrand = range.getTargetStrand() == TargetStrand.POSITIVE_STRAND;

    if ((queryRange.getEnd() - queryRange.getStart())
        != (targetRange.getEnd() - targetRange.getStart())) {
      return RegionType.ALIGNMENT_REQUIRED;
    }

    String querySeq = queryFasta.get(queryRange.getReferenceName(),
        queryRange.getStart(), queryRange.getEnd());

    String targetSeq = targetFasta.get(targetRange.getReferenceName(),
        targetRange.getStart(), targetRange.getEnd());

    if (Fasta.MISSING_CHROMOSOME.equals(targetSeq)) {
      return RegionType.UNKNOWN_REGION_TYPE;
    }

    if (!GenomeWarpUtils.isValidDna(querySeq) || !GenomeWarpUtils.isValidDna(targetSeq)) {
      // This type will be filtered out
      return RegionType.UNKNOWN_REGION_TYPE;
    }

    if ((isPositiveStrand && querySeq.equals(targetSeq))
        || (!isPositiveStrand && querySeq.equals(SequenceUtil.reverseComplement(targetSeq)))) {
      return RegionType.IDENTICAL;
    }

    return RegionType.MISMATCHED_BASES;
  }


  /**
   * Given a series of range's transformations, classifies each range into one of 5 types.
   *
   * <p> This function takes a homologous range, which dictates the query range and
   * the target range resulting from a lift over, then classifies it as ALIGNMENT_REQUIRED,
   * IDENTICAL, or MISMATCHED_BASES. UNKNOWN_REGION_TYPES are removed, as they indicate
   * non-DNA bases.
   *
   * @param queryFasta the reference query genome
   * @param targetFasta the reference target genome
   * @param ranges a list of HomologousRanges encoding transformations
   * @return a List of HomologousRanges, including region type
   */
  public static List<HomologousRange> analyzeRegions(Fasta queryFasta, Fasta targetFasta,
      List<HomologousRange> ranges) {
    for (ListIterator<HomologousRange> iterator = ranges.listIterator(); iterator.hasNext();) {
      HomologousRange currRange = iterator.next();
      RegionType regionType = getRegionType(currRange, queryFasta, targetFasta);

      if (regionType == RegionType.UNKNOWN_REGION_TYPE) {
        iterator.remove();
      } else {
        iterator.set(HomologousRange.newBuilder(currRange).setRegionType(regionType).build());
      }
    }

    return ranges;
  }

  private static SortedMap<String, List<GenomeRange>> performLiftOver(List<GenomeRange> inRange) {

    SortedMap<String, List<GenomeRange>> toReturn = new TreeMap<>();

    LiftOver liftOverTool = new LiftOver(new File(ARGS.liftOverChainPath));

    int all = inRange.size();
    int current = 0;
    int currentMarker = 1;

    for (GenomeRange currRange : inRange) {
      int oneBasedClosedStart = (int) currRange.getStart() + 1;
      int oneBasedClosedEnd = (int) currRange.getEnd();
      if ((long) oneBasedClosedStart != currRange.getStart() + 1
          || (long) oneBasedClosedEnd != currRange.getEnd()) {
        throw new IllegalArgumentException("performLiftOver: cannot specify ranges which are "
            + "greater than Integer.MAX_VALUE");
      }

      // These intervals are 1-based, close ended
      final Interval currInterval = new Interval(currRange.getChromosome(), oneBasedClosedStart,
          oneBasedClosedEnd, !currRange.isPositiveStrand(), currRange.getName());

      Interval liftedInterval = liftOverTool.liftOver(currInterval, ARGS.minMatch);
      if (liftedInterval == null) {
        logger.log(Level.FINE, "failed to liftover an interval");
      } else {
        // Change from one based to 0 based
        GenomeRange toAdd = new GenomeRange(liftedInterval.getSequence(),
            liftedInterval.getStart() - 1, liftedInterval.getEnd(), liftedInterval.getName(),
            liftedInterval.isPositiveStrand());
        if (!toReturn.containsKey(liftedInterval.getSequence())) {
          toReturn.put(liftedInterval.getSequence(), new ArrayList<GenomeRange>());
        }
        toReturn.get(liftedInterval.getSequence()).add(toAdd);
      }

      if ((100 * current++) / all > currentMarker) {
        logger.log(Level.INFO, String.format("performing liftover: at %d%%", currentMarker++));
      }
    }

    return toReturn;
  }

  private static void writeHeader(List<String> header, PrintWriter writer) {
    for (String line : header) {
      writer.println(line);
    }
  }

  private static Map<Integer, List<HomologousRange>> saveToFile(
      Map<HomologousRange, List<String>> groupedVariants, List<String> header,
      Map<Integer, List<String>> queryChr, Map<Integer, List<String>> targetChr)
      throws FileNotFoundException {
    int parts = 0;
    int headerSize = header.size();
    int total = headerSize;
    Map<Integer, List<HomologousRange>> toReturn = new HashMap<>();
    PrintWriter vcfWriter = new PrintWriter(ARGS.workDir + VCF_PART_NAME + (++parts));
    writeHeader(header, vcfWriter);

    List<HomologousRange> currRegionList = new ArrayList<>();
    List<HomologousRange> keyList = new ArrayList<>(groupedVariants.keySet());

    // Store which chromosomes to preload for each part
    List<String> currQueryChrList = new ArrayList<>();
    List<String> currTargetChrList = new ArrayList<>();

    Collections.sort(keyList, new GenomeWarpUtils.HomologousRangeComparator());
    for (HomologousRange currRegion : keyList) {
      List<String> currList = groupedVariants.get(currRegion);

      total += currList.size();
      // If we have too many lines, add to next file. However, never
      // create files with only a header
      if (total >= ARGS.vcfChunkSize && total != headerSize) {
        // Save which regions are in this vcf chunk
        toReturn.put(Integer.valueOf(parts), currRegionList);
        currRegionList = new ArrayList<>();

        // Save the chromsomes seen for this part
        queryChr.put(Integer.valueOf(parts), currQueryChrList);
        targetChr.put(Integer.valueOf(parts), currTargetChrList);
        currQueryChrList = new ArrayList<>();
        currTargetChrList = new ArrayList<>();

        // Close old writer and open new one
        vcfWriter.close();
        vcfWriter = new PrintWriter(ARGS.workDir + VCF_PART_NAME + (++parts));
        writeHeader(header, vcfWriter);
        total = currList.size() + headerSize;
      }

      String currQueryChr = currRegion.getQueryRange().getReferenceName();
      String currTargetChr = currRegion.getTargetRange().getReferenceName();
      if (!currQueryChrList.contains(currQueryChr)) {
        currQueryChrList.add(currQueryChr);
      }
      if (!currTargetChrList.contains(currTargetChr)) {
        currTargetChrList.add(currTargetChr);
      }

      currRegionList.add(currRegion);
      for (String currLine : currList) {
        vcfWriter.println(currLine);
      }
    }
    // Fence post
    queryChr.put(Integer.valueOf(parts), currQueryChrList);
    targetChr.put(Integer.valueOf(parts), currTargetChrList);
    toReturn.put(Integer.valueOf(parts), currRegionList);
    vcfWriter.close();

    return toReturn;
  }

  private static Map<String, String> createContigEntry(String id, long length, String assembly,
      String species) {
    // We use linked hash map to maintain insertion order
    Map<String, String> toReturn = new LinkedHashMap<>();
    toReturn.put("ID", id);
    toReturn.put("length", String.valueOf(length));
    toReturn.put("assembly", assembly);
    toReturn.put("species", species);

    return toReturn;
  }

  private static VCFHeader warpHeader(VCFHeader in, Map<String, Long> namesAndLength)
      throws IllegalArgumentException {
    Set<VCFHeaderLine> newLines = new HashSet<>();
    boolean hasSource = false;
    VCFHeaderVersion version = DEFAULT_VCF_VERSION;

    for (VCFHeaderLine line : in.getMetaDataInInputOrder()) {
      if (line.getKey().equals("reference")) {
        newLines.add(new VCFHeaderLine(line.getKey(), ARGS.refTargetFASTA));
      } else if (line.getKey().equals("source")) {
        newLines.add(new VCFHeaderLine(line.getKey(), line.getValue() + "_and_"
            + GENOME_WARP_VERSION));
        hasSource = true;
      } else if (line.getKey().equals("fileformat")) {
        version = VCFHeaderVersion.toHeaderVersion(line.getValue());
        if (version == null) {
          throw new IllegalArgumentException("malformed version: " + line.getValue());
        }
      } else if (line.getKey().equals(VCFConstants.CONTIG_HEADER_KEY)) {
        continue;
      } else {
        newLines.add(line);
      }
    }

    if (!hasSource) {
      newLines.add(new VCFHeaderLine("source", GENOME_WARP_VERSION));
    }

    // Add contigs
    int i = 0;
    for (Map.Entry<String, Long> entry : namesAndLength.entrySet()) {
      String currName = entry.getKey();
      long chrSize = entry.getValue();

      newLines.add(new VCFContigHeaderLine(VCFHeaderLine.toStringEncoding(
          createContigEntry(currName, chrSize, ARGS.targetAssembly, ARGS.species)),
          version, VCFConstants.CONTIG_HEADER_KEY, i++));
    }

    return new VCFHeader(newLines, in.getSampleNamesInOrder());
  }

  /**
   * Reads in the initial BED file, prepares it for liftover,
   * then performs liftover on BED, then classifies regions
   * based on pre- and post- liftover regions.
   */
  private static List<HomologousRange> processAndSaveBed(String inputBed, String inputVcf,
      Fasta queryFasta,
      Fasta targetFasta) {
    // Open necessary readers
    BufferedReader bedReader = null;
    boolean simplifiedPreprocessing = ARGS.simplifiedRegionsPreprocessing;
    try {
      logger.log(Level.INFO, "Reading BED");
      bedReader = Files.newBufferedReader(Paths.get(inputBed), UTF_8);
    } catch (IOException ex) {
      GenomeWarpUtils.fail(logger, "failed to parse input BED/FASTA/VCF file(s): "
          + ex.getMessage());
    }

    // Takes the input BED and splits at non-DNA characters
    logger.log(Level.INFO, "Split DNA at non-DNA characters");
    SortedMap<String, List<GenomeRange>> dnaOnlyInputBEDPerChromosome = null;
    try {
      if ((dnaOnlyInputBEDPerChromosome = GenomeRangeUtils.splitAtNonDNA(queryFasta, bedReader)) == null) {
        GenomeWarpUtils.fail(logger, "failed to generate reference genome");
      }
    } catch (IOException ex) {
      GenomeWarpUtils.fail(logger, "failed to read from input bed: " + ex.getMessage());
    }

    /**
     * PRE PROCESSING
     */
    List<GenomeRange> queryBED = (simplifiedPreprocessing) ?
        GenomeRangeUtils.generateQueryBEDWithSimplifiedPreprocessing(dnaOnlyInputBEDPerChromosome):
        GenomeRangeUtils
            .generateQueryBEDWithImprovedPreprocessing(dnaOnlyInputBEDPerChromosome, inputVcf,
                ARGS.bedWindowSize);

    /**
     * LIFTOVER
     **/
    logger.log(Level.INFO, String.format("Performing liftover on %d ranges", queryBED.size()));
    SortedMap<String, List<GenomeRange>> liftedBEDPerChromosome = performLiftOver(queryBED);
    List<GenomeRange> targetBED = new ArrayList<>();

    /**
     * POST LIFTOVER
     **/
    for (List<GenomeRange> liftedBEDChr: liftedBEDPerChromosome.values()) {
      Collections.sort(liftedBEDChr);

      logger.log(Level.INFO, "Removing overlap");

      // Remove overlap in processed BED file
      targetBED.addAll(GenomeRangeUtils.omitOverlap(liftedBEDChr));
    }

    logger.log(Level.INFO, "Starting join regions and analysis");

    // Create joined structure for query/target BED files
    List<HomologousRange> joinedRegions = GenomeRangeUtils.joinRegions(queryBED, targetBED);
    List<HomologousRange> namedRegions = analyzeRegions(queryFasta, targetFasta, joinedRegions);

    if (ARGS.intermediateFiles) {
      try {
        PrintWriter out = new PrintWriter(Files.newBufferedWriter(Paths.get(ARGS.workDir
            + "/out.bed.annotated"), UTF_8));
        for (HomologousRange currRegion : namedRegions) {
          out.println(GenomeWarpUtils.homologousToString(currRegion));
        }
        out.close();
      } catch (IOException ex) {
        GenomeWarpUtils.fail(logger, "failed to write out annotated regions: "
            + ex.getMessage());
      }
    }

    try {
      bedReader.close();
    } catch (IOException ex) {
      GenomeWarpUtils.fail(logger, "failed to close opened buffers: " + ex.getMessage());
    }

    return namedRegions;
  }

  public static void main(String[] args) {
    JCommander parser = new JCommander(ARGS);
    try {
      parser.parse(args);
    } catch (ParameterException e) {
      parser.usage();
      GenomeWarpUtils.fail(logger, "failed to parse command line args: " + e.getMessage());
    }

    List<String> headerStrings = null;
    List<HomologousRange> namedRegions = null;
    Fasta queryFasta = null, targetFasta = null;

    // Solve Issue #2 - Check if VCS is a gVCF and if so extract 1) variant-only VCF, 2) BED file
    boolean haveGvcf = ARGS.rawQueryGvcf != null;
    final String queryVcfToProcess;
    final String queryBedToProcess;
    // Query input flag validation.
    if (haveGvcf && (ARGS.rawQueryVcf != null || ARGS.rawQueryBed != null)) {
      GenomeWarpUtils.fail(logger,
          "Arguments (--raw_query_vcf, --raw_query_bed) and --raw_query_gvcf are mutually exclusive");
    }
    if ((ARGS.rawQueryVcf == null) != (ARGS.rawQueryBed == null)) {
      GenomeWarpUtils.fail(logger,
          "Either both or neither of --raw_query_vcf and --raw_query_bed must be specified");
    }
    if (ARGS.rawQueryGvcf == null && ARGS.rawQueryVcf == null) {
      GenomeWarpUtils.fail(logger,
          "Either (--raw_query_vcf, --raw_query_bed) or --raw_query_gvcf must be specified");
    }
    if (haveGvcf) {
      queryVcfToProcess = ARGS.workDir + File.separator + "from_gvcf.vcf";
      queryBedToProcess = ARGS.workDir + File.separator + "from_gvcf.bed";
      logger.log(Level.INFO, "Checking and processing gVCF");
      if (!GvcfToVcfAndBed.saveVcfAndBedFromGvcf(ARGS.rawQueryGvcf, queryVcfToProcess,
          queryBedToProcess)) {
        GenomeWarpUtils.fail(logger, "Failed to read gVCF/write VCF or BED files");
      }
      logger.log(Level.INFO, "Recognized gVCF format - using extracted VCF and BED files as input");
    } else {
      // regular VCF
      queryVcfToProcess = ARGS.rawQueryVcf;
      queryBedToProcess = ARGS.rawQueryBed;
    }      
 
    logger.log(Level.INFO, "Creating FASTA structure and jump table");
    try {
      headerStrings = retrieveVcfHeader(queryVcfToProcess);
      queryFasta = new Fasta(ARGS.refQueryFASTA);
      targetFasta = new Fasta(ARGS.refTargetFASTA);
    } catch (IOException ex) {
      GenomeWarpUtils.fail(logger, "failed initial setup: " + ex.getMessage());
    }

    if (!ARGS.onlyGenomeWarp) {
      namedRegions = processAndSaveBed(queryBedToProcess, queryVcfToProcess, queryFasta, targetFasta);
    } else {
      logger.log(Level.INFO, "Region data from files");
      BufferedReader regionsFile = null;
      try {

        regionsFile = Files.newBufferedReader(Paths.get(ARGS.workDir
            + "/out.bed.annotated"), UTF_8);

        // Load Regions from file
        String line;
        namedRegions = new ArrayList<>();
        while ((line = regionsFile.readLine()) != null) {
          namedRegions.add(GenomeWarpUtils.homologousFromString(line));
        }

        regionsFile.close();
      } catch (IOException ex) {
        GenomeWarpUtils.fail(logger, "failed to parse input FASTA/Region file(s): "
            + ex.getMessage());
      }
    }

    // Read variants into variant contexts and group into regions
    logger.log(Level.INFO, "Finished preloading, grouping variants into regions");
    Map<HomologousRange, List<String>> groupedVariants = null;
    Map<Integer, List<HomologousRange>> mapping = null;
    Map<Integer, List<String>> queryChr = new HashMap<>();
    Map<Integer, List<String>> targetChr = new HashMap<>();
    try {
      BufferedReader vcfFile =
          Files.newBufferedReader(Paths.get(queryVcfToProcess), UTF_8);
      groupedVariants = GenomeWarpUtils.associateVariantsWithRange(namedRegions, vcfFile);
      logger.log(Level.INFO, "Saving into discrete files of ~50000 variants");
      mapping = saveToFile(groupedVariants, headerStrings, queryChr, targetChr);
    } catch (IOException ex) {
      GenomeWarpUtils.fail(logger, "Failed to open vcf file");
    }

    /**
     * Begin actual GenomeWarping
     */
    VCFFileReader vcfReader =
        new VCFFileReader(new File(queryVcfToProcess), false);
    VCFHeader vcfHeader = vcfReader.getFileHeader();
    ArrayList<String> vcfSampleNames = vcfHeader.getSampleNamesInOrder();
    if (vcfSampleNames.size() == 0) {
      GenomeWarpUtils.fail(logger, "input VCF file has no callset groups");
    }
    final boolean keepRefRefCalls = ARGS.keepHomozygousReferenceCalls;

    // Minimize memory usage by operating over 50 kilobase chunks
    PrintWriter outBed = null;
    FileOutputStream outVcf = null;
    try {
      outBed = new PrintWriter(Files.newBufferedWriter(Paths.get(ARGS.outputRegionsFile),
          UTF_8));
      outVcf = new FileOutputStream(new File(ARGS.outputVariantsFile));
    } catch (IOException ex) {
      GenomeWarpUtils.fail(logger, "failed to write variants and confidently-called regions to file: "
          + ex.getMessage());
    }

    logger.log(Level.INFO, "Warping VCF header");
    try {
      vcfHeader = warpHeader(vcfHeader, targetFasta.getReferenceNameAndLengthInInputOrder());
    } catch (IllegalArgumentException ex) {
      GenomeWarpUtils.fail(logger, "warpHeader failure: " + ex.getMessage());
    }

    int parts = mapping.size();
    List<GenomeRange> finalRanges = new ArrayList<>();
    for (int part = 1; part <= parts; part++) {
      List<Variant> finalVariants = new ArrayList<>();
      // Display progress logs
      logger.log(Level.INFO, String.format("Processing part %d of %d", part, parts));

      // Convert to variants
      List<HomologousRange> currRegions = mapping.get(Integer.valueOf(part));
      List<Variant> variants = null;
      try {
        File thisFile = new File(ARGS.workDir + VCF_PART_NAME + part);
        thisFile.deleteOnExit();
        variants = VcfToVariant.convertVcfToVariant(thisFile, false);
      } catch (IOException ex) {
        GenomeWarpUtils.fail(logger, "Failed to open partial vcf");
      }

      if (ARGS.removeInfoField) {
        variants = GenomeWarpUtils.removeInfoField(variants);
      }

      // Create mappings from range to variants
      Map<HomologousRange, List<Variant>> intermediate = GenomeWarpUtils.associateVariantsWithRange(
          currRegions, variants);

      // Preload FASTA
      try {
        queryFasta.preload(queryChr.get(Integer.valueOf(part)));
        targetFasta.preload(targetChr.get(Integer.valueOf(part)));
      } catch (IOException ex) {
        GenomeWarpUtils.fail(logger, String.format("Failed to preload FASTA: %s", ex.getCause()));
      }

      List<HomologousRange> keyList = new ArrayList<>(intermediate.keySet());
      Collections.sort(keyList, new GenomeWarpUtils.HomologousRangeComparator());
      for (HomologousRange currRegion : keyList) {
        List<Variant> currVariants = intermediate.get(currRegion);
        List<Variant> targetVariants = TransformVariantsInRange.transformQueryToTargetVariants(
            currRegion, currVariants, vcfSampleNames, queryFasta, targetFasta);

        if (targetVariants == null) {
          continue;
        }

        /**
         * Add BED to be printed
         */
        String regionChr = currRegion.getTargetRange().getReferenceName();
        Long regionStart = currRegion.getTargetRange().getStart();
        Long regionEnd = currRegion.getTargetRange().getEnd();
        finalRanges.add(new GenomeRange(regionChr, regionStart, regionEnd));

        // Save variants
        for (Variant currVariant : targetVariants) {
          if (GenomeWarpUtils.hasVariation(currVariant) || keepRefRefCalls) {
            finalVariants.add(currVariant);
          }
        }
      }

      // Sort and print variants for this part
      Comparator<Variant> comp = new GenomeWarpUtils.VariantComparator();
      if (finalVariants.size() > 0) {
        Collections.sort(finalVariants, comp);
        VariantToVcf.convertVariantToVcf(vcfHeader, finalVariants, outVcf, part == 1);
      }
    }

    logger.log(Level.INFO, "Sorting final BEDs in preparation for output");
    Collections.sort(finalRanges);
    finalRanges = GenomeRangeUtils.mergeOverlaps(finalRanges);
    for (GenomeRange currRange : finalRanges) {
      outBed.println(String.format("%s\t%s\t%s", currRange.getChromosome(),
          Long.toString(currRange.getStart()), Long.toString(currRange.getEnd())));
    }

    /**
     * Finished all file reading, close opened buffers
     **/
    try {
      outBed.close();
      outVcf.close();
    } catch (IOException ex) {
      GenomeWarpUtils.fail(logger, "failed to close opened buffers: " + ex.getMessage());
    }

  }
}
