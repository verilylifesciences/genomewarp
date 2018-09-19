package com.verily.genomewarp.utils;

import com.google.genomics.v1.Range;
import com.verily.genomewarp.GenomeWarpSerial;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.TargetStrand;
import com.verily.genomewarp.utils.GenomeWarpUtils.HomologousRangeComparator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.SortedMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;

public class GenomeRangeUtils {

  private static final Logger logger = Logger.getLogger(GenomeRangeUtils.class.getName());

  /**
   * Splits a given set of genome ranges at non-DNA characters
   *
   * <p> This function reads in ranges using a BufferedReader and splits
   * the given range along non-strict DNA characters (non ACTGactg).
   *
   * @param refFasta the reference genome, containing the string bases
   * @param bed a BufferedReader containing the range reads
   * @return a hashmap of chromosome to arraylist of GenomeRanges only containing valid DNA characters
   */
  public static SortedMap<String, List<GenomeRange>> splitAtNonDNA(Fasta refFasta, BufferedReader bed)
      throws IOException {

    String line;
    String pastChr = "";
    SortedMap<String, List<GenomeRange>> matchesPerChr = new TreeMap<>();
    List<GenomeRange> currMatches = null;
    while ((line = bed.readLine()) != null) {
      String[] range = line.trim().split("\\s+");
      if (range.length < 3) {
        GenomeWarpUtils.fail(logger, "input BED file has less than 3 columns");
      }

      String chr = range[0];
      long start = Long.parseLong(range[1]);
      long end = Long.parseLong(range[2]);

      if (!chr.equals(pastChr)) {
        logger.log(Level.INFO, String.format("Now processing %s", chr));
        pastChr = chr;
        currMatches = new ArrayList<>();
        matchesPerChr.put(chr, currMatches);
      }

      String chrSubSeq = refFasta.get(chr, start, end);
      Matcher matchDNASeq = GenomeWarpUtils.dnaSequence.matcher(chrSubSeq);

      // Use the regex matcher to add subsequences which are
      // valid DNA base runs to a match list
      int splitNum = 0;
      while (matchDNASeq.find()) {
        GenomeRange currRange = new GenomeRange(chr, start + matchDNASeq.start(),
            start + matchDNASeq.end());
        currMatches.add(currRange);
        splitNum++;
      }
      if (splitNum > 1) {
        logger.log(Level.INFO,
          String.format("%s (%d, %d) was split into %d regions", chr, start, end, splitNum));
      }
    }

    SortedMap<String, List<GenomeRange>> nonEmptyMatchesPerChr = new TreeMap<>();
    for (SortedMap.Entry<String, List<GenomeRange>> entry: matchesPerChr.entrySet())
    {
      if (!entry.getValue().isEmpty()) {
        nonEmptyMatchesPerChr.put(entry.getKey(), entry.getValue());
      }
    }

    return nonEmptyMatchesPerChr;
  }

  /**
   * Generates set of genome ranges per chromosome from the variants in the VCF file
   *
   * For each variant in the input VCF, it will create a region corresponding to the this variant.
   *
   * @param vcfReader the input VCF file reader
   * @return a map of chromosome to the list of GenomeRanges corresponding to the variants in the input VCF
   */
  public static SortedMap<String, List<GenomeRange>> generateBEDFromVCF(VCFFileReader vcfReader) {
    SortedMap<String, List<GenomeRange>> toReturn = new TreeMap<>();
    for (final VariantContext var: vcfReader) {
      // TODO: check that start and end positions are valid
      // htsjdk start/end are both 1-based closed
      GenomeRange curr = new GenomeRange(var.getChr(),
              var.getStart() - 1, var.getEnd());
      if (!toReturn.containsKey(curr.getChromosome())) {
        toReturn.put(curr.getChromosome(), new ArrayList<GenomeRange>());
      }
      toReturn.get(curr.getChromosome()).add(curr);
    }
    return toReturn;
  }

  /**
   * Splits the given region and returns a list subregions of at most specified length.
   *
   * @param region GenomeRange to split
   * @param windowSize Maximum length of the regions in result
   * @return list of subregions generated from the input region
   */
  public static List<GenomeRange> splitRegion(GenomeRange region, int windowSize) {
    List<GenomeRange> toReturn = new ArrayList<>();
    long pos = region.getStart();
    for (; pos + windowSize <= region.getEnd(); pos += windowSize) {
      toReturn.add(new GenomeRange(region.getChromosome(), pos, pos + windowSize));
    }
    if (pos != region.getEnd()) {
      toReturn.add(new GenomeRange(region.getChromosome(), pos, region.getEnd()));
    }
    return toReturn;
  }

  /**
   * Filters out variant regions not covered by confident regions.
   *
   * Only variant regions completely covered by one of confident regions will be included in the
   * output.
   *
   * Inputs are expected to include only genome ranges from the same chromosome.
   *
   * @param queryBED sorted list of GenomeRanges from the query BED file coming from the same chromosome
   * @param fromVcfBED sorted list of GenomeRanges generated from VCF file coming from the same chromosome
   * @return filtered sorted list of GenomeRanges from fromVcfBED
   */
  public static List<GenomeRange> filterOutNotCoveredVariants(
      List<GenomeRange> queryBED, List<GenomeRange> fromVcfBED) {
    List<GenomeRange> toReturn = new ArrayList<>();
    ListIterator<GenomeRange> queryIt = queryBED.listIterator();
    ListIterator<GenomeRange> vcfIt = fromVcfBED.listIterator();
    GenomeRange queryRegion = GenomeWarpUtils.nextOrNull(queryIt);
    GenomeRange vcfRegion = GenomeWarpUtils.nextOrNull(vcfIt);

    while (queryRegion != null && vcfRegion != null) {
      if (!queryRegion.getChromosome().equals(vcfRegion.getChromosome())) {
        GenomeWarpUtils.fail(logger, "inputs must contain the data from a single chromosome");
      }

      if (queryRegion.getEnd() <= vcfRegion.getStart()) {
        queryRegion = GenomeWarpUtils.nextOrNull(queryIt);
      } else if (vcfRegion.getEnd() <= queryRegion.getStart()) {
        vcfRegion = GenomeWarpUtils.nextOrNull(vcfIt);
      } else {
        if (queryRegion.includes(vcfRegion)) {
          toReturn.add(vcfRegion);
        }
        vcfRegion = GenomeWarpUtils.nextOrNull(vcfIt);
      }
    }
    return toReturn;
  }

  /**
   * Pads ranges in the input by VARIANT_CONTEXT_SIZE bps from each side.
   *
   * @param ranges list of GenomeRanges
   * @return list of padded GenomeRanges
   */
  public static List<GenomeRange> generatePaddedGenomeRanges(
      List<GenomeRange> ranges) {
    List<GenomeRange> toReturn = new ArrayList<>();
    for (GenomeRange range: ranges) {
      toReturn.add(new GenomeRange(range.getChromosome(),
          range.getStart() - GenomeWarpUtils.VARIANT_CONTEXT_SIZE,
          range.getEnd() + GenomeWarpUtils.VARIANT_CONTEXT_SIZE));
    }
    return toReturn;
  }

  /**
   * Generates set of genome ranges by merging regions from input BED and regions generated from VCF variants
   *
   * The merging strategy is the following: any VCF region that doesn't overlap query BED will omitted.
   * Any query region that overlap a variant are truncated to the start of the variant region.
   * Variant regions that overlap any query region will be included in the result.
   * All regions from query bed are splitted to have at most bedWindowSize length.
   *
   * Inputs are expected to include only genome ranges from the same chromosome.
   *
   * @param queryBED sorted list of GenomeRanges from the query BED file coming from the same chromosome
   * @param fromVcfBED sorted list of GenomeRanges generated from VCF file coming from the same chromosome
   * @param windowSize regions from queryBED will be splitted to subregions of at most this size
   * @return sorted list of genome ranges merged from queryBED and fromVcfBED
   */
  public static List<GenomeRange> mergeRegionsFromQueryBEDAndVariants(
      List<GenomeRange> queryBED, List<GenomeRange> fromVcfBED, int windowSize) {
    List<GenomeRange> mergedBed = new ArrayList<>();
    ListIterator<GenomeRange> queryIt = queryBED.listIterator();
    ListIterator<GenomeRange> vcfIt = fromVcfBED.listIterator();
    GenomeRange queryRegion = GenomeWarpUtils.nextOrNull(queryIt);
    GenomeRange vcfRegion = GenomeWarpUtils.nextOrNull(vcfIt);

    while (queryRegion != null) {
      if (vcfRegion == null) {
          mergedBed.addAll(splitRegion(queryRegion, windowSize));
          queryRegion = GenomeWarpUtils.nextOrNull(queryIt);
          continue;
      }

      if (!queryRegion.getChromosome().equals(vcfRegion.getChromosome())) {
        GenomeWarpUtils.fail(logger,"inputs must contain data from a single chromosome");
      }

      if (queryRegion.getEnd() <= vcfRegion.getStart()) {
        mergedBed.addAll(splitRegion(queryRegion, windowSize));
        queryRegion = GenomeWarpUtils.nextOrNull(queryIt);
      }
      else if (vcfRegion.getEnd() <= queryRegion.getStart()) {
        vcfRegion = GenomeWarpUtils.nextOrNull(vcfIt);
      }
      else {
        if (queryRegion.getStart() < vcfRegion.getStart()) {
          mergedBed.addAll(splitRegion(
              new GenomeRange(queryRegion.getChromosome(), queryRegion.getStart(),
                  vcfRegion.getStart()), windowSize
          ));
        }
        mergedBed.add(queryRegion.getIntersection(vcfRegion));
        if (queryRegion.getEnd() > vcfRegion.getEnd()) {
          queryRegion = new GenomeRange(queryRegion.getChromosome(), vcfRegion.getEnd(),
              queryRegion.getEnd());
          vcfRegion = GenomeWarpUtils.nextOrNull(vcfIt);
        } else {
          queryRegion = GenomeWarpUtils.nextOrNull(queryIt);
        }
      }
    }
    return mergedBed;
  }

  /**
   * Merges overlapping region (including consecutive regions without a gap) together
   *
   * @param toMerge sorted list of GenomeRanges to be merged
   * @return sorted list of non-overlapping GenomeRanges
   */
  public static List<GenomeRange> mergeOverlaps(List<GenomeRange> toMerge) {
    GenomeRange maxRange = null;
    List<GenomeRange> toReturn = new ArrayList<>();
    for (GenomeRange currRange: toMerge) {
      if (maxRange != null) {
        if (currRange.getChromosome().equals(maxRange.getChromosome())
            && currRange.compareTo(maxRange) < 0) {
          GenomeWarpUtils.fail(logger, "input regions are not sorted by position");
        }

        if (!maxRange.getChromosome().equals(currRange.getChromosome())
            || maxRange.getEnd() < currRange.getStart()) {
          toReturn.add(maxRange);
          maxRange = currRange;
        } else if (maxRange.getEnd() < currRange.getEnd()) {
          maxRange = new GenomeRange(maxRange.getChromosome(), maxRange.getStart(), currRange.getEnd());
        }
      } else {
        maxRange = currRange;
      }
    }
    if (maxRange != null) {
      toReturn.add(maxRange);
    }
    return toReturn;
  }

  /**
   * Applies final preprocessing operations to the GenomeRanges list
   *
   * @param inBED sorted list of GenomeRanges
   * @return list of GenomeRanges
   */
  public static List<GenomeRange> massageBED(List<GenomeRange> inBED) {
    List<GenomeRange> toReturn = new ArrayList<>();

    int i = 1;
    for (GenomeRange currRange : inBED) {
      String lineName = currRange.getChromosome() + "." + Integer.toString(i);
      currRange.setName(lineName);

      toReturn.add(currRange);
      i++;
    }

    return toReturn;
  }

  /**
   * Removes all BED regions that overlap with others
   *
   * <p> This function iterates through all ArrayLists, while maintaining the
   * maximum values of end attained, ensuring no future BEDS intersect this
   * range. Furthermore, this function requires the input BED be sorted
   * by chromosome and starting position.
   *
   * Inputs are expected to include only genome ranges from the same chromosome.
   *
   * @param inBED the sorted array of BED ranges to be analyzed. All ranges must come from the same chromosome.
   * @return an arraylist omitting any BED regions which overlap
   */
  public static List<GenomeRange> omitOverlap(List<GenomeRange> inBED) {
    GenomeRange prevRange = null;
    GenomeRange maxRange = null;
    boolean omitPrev = true;

    List<GenomeRange> nonOverlappingBED = new ArrayList<>();

    for (GenomeRange currRange : inBED) {

      // Check that input BED is sorted
      if (prevRange != null) {
        if (!currRange.getChromosome().equals(prevRange.getChromosome())) {
          GenomeWarpUtils.fail(logger,"found ranges from different chromosomes");
        }
        if (currRange.getStart() < prevRange.getStart()) {
          GenomeWarpUtils.fail(logger, "output BED of liftover is not sorted by position");
        }
      }

      boolean omitCurr = currRange.overlaps(maxRange);
      omitPrev = omitPrev || omitCurr;

      if (!omitPrev) {
        nonOverlappingBED.add(prevRange);
      }

      if (prevRange != null) {
        if (maxRange == null) {
          GenomeWarpUtils.fail(logger, "something is amiss in omitOverlap");
        }
        if (maxRange.getEnd() < currRange.getEnd()) {
          maxRange = new GenomeRange(maxRange.getChromosome(), maxRange.getStart(),
              currRange.getEnd());
        }
      } else {
        maxRange = currRange;
      }

      omitPrev = omitCurr;
      prevRange = currRange;
    }

    // Handle fencepost
    if (!omitPrev) {
      nonOverlappingBED.add(prevRange);
    }

    return nonOverlappingBED;
  }

  /**
   * Given an array of query and target regions we join them by name.
   *
   * @param queryBED an List of query BED ranges
   * @param targetBED an List of target BED ranges
   * @return an List of Homologous range protos
   */
  public static List<HomologousRange> joinRegions(List<GenomeRange> queryBED,
      List<GenomeRange> targetBED) {
    List<HomologousRange> associate = new ArrayList<>();
    Map<String, GenomeRange> mappedQueries = new HashMap<>();
    Set<String> seenTargetNames = new HashSet<String>();

    // Add all query ranges. Fail if we encounter multiple query BEDS
    // with the same name
    for (GenomeRange currRange : queryBED) {
      if (mappedQueries.containsKey(currRange.getName())) {
        GenomeWarpUtils.fail(logger, "found duplicated BED names in query BED");
      }

      mappedQueries.put(currRange.getName(), currRange);
    }

    // Add all target ranges. Fail if we encounter multiple target BEDS
    // with the same name
    for (GenomeRange currTarget : targetBED) {
      if (seenTargetNames.contains(currTarget.getName())) {
        GenomeWarpUtils.fail(logger, "found duplicated BED names in target BED");
      }
      seenTargetNames.add(currTarget.getName());

      GenomeRange query = mappedQueries.get(currTarget.getName());
      if (query != null) {
        associate.add(HomologousRange.newBuilder().setQueryRange(
            Range.newBuilder()
            .setReferenceName(query.getChromosome())
            .setStart(query.getStart())
            .setEnd(query.getEnd()))
          .setTargetRange(
            Range.newBuilder()
            .setReferenceName(currTarget.getChromosome())
            .setStart(currTarget.getStart())
            .setEnd(currTarget.getEnd()))
          .setTargetStrand(currTarget.isPositiveStrand()
            ? TargetStrand.POSITIVE_STRAND : TargetStrand.NEGATIVE_STRAND).build());
      }
    }

    Collections.sort(associate, new HomologousRangeComparator());
    return associate;
  }

  /**
   * Performs preprocessing of the input confident regions.
   *
   * @param inputBEDPerChromosome map of chromosome to the list of confident regions. Confident regions must not contain non-DNA characters.
   * @return list of the GenomeRanges ready for lift over to the target genome
   */
  public static List<GenomeRange> generateQueryBEDWithSimplifiedPreprocessing(
      SortedMap<String, List<GenomeRange>> inputBEDPerChromosome) {
    List<GenomeRange> queryBED = new ArrayList<>();
    for (String chromosome: inputBEDPerChromosome.keySet()) {
      List<GenomeRange> inputBEDChr = inputBEDPerChromosome.get(chromosome);
      Collections.sort(inputBEDChr);

      logger.log(Level.INFO, String.format("Massaging %d BED record(s) from chromosome %s",
          inputBEDChr.size(), chromosome));
      queryBED.addAll(massageBED(inputBEDChr));
    }
    return queryBED;
  }

  /**
   * Performs improved preprocessing of the input confident regions
   *
   * Splits the confident regions to do local remapping around variants in the input VCF to maximize
   * recall of variants.
   *
   * @param inputBEDPerChromosome map of chromosome to the list of confident regions. Confident regions must not contain non-DNA characters.
   * @param inputVcf input VCF file
   * @param bedWindowSize regions from intput BED will be splitted to subregions of at most this size
   * @return list of the GenomeRanges ready for lift over to the target genome
   */
  public static List<GenomeRange> generateQueryBEDWithImprovedPreprocessing(
      SortedMap<String, List<GenomeRange>> inputBEDPerChromosome, String inputVcf,
      int bedWindowSize) {
    List<GenomeRange> queryBED = new ArrayList<>();
    VCFFileReader vcfReader = new VCFFileReader(new File(inputVcf), false);

    logger.log(Level.INFO, "Generating regions from variants");
    SortedMap<String, List<GenomeRange>> fromVcfBEDPerChromosome = generateBEDFromVCF(vcfReader);

    for (String chromosome: inputBEDPerChromosome.keySet()) {
      List<GenomeRange> inputBEDChr = inputBEDPerChromosome.get(chromosome);
      Collections.sort(inputBEDChr);

      List<GenomeRange> fromVcfBEDChr = fromVcfBEDPerChromosome.get(chromosome);
      if (fromVcfBEDChr == null) {
        continue;
      }
      Collections.sort(fromVcfBEDChr);

      logger.log(Level.INFO,
          "Filtering out VCF regions not covered by confident regions");
      List<GenomeRange> filteredVcfBEDChr = filterOutNotCoveredVariants(inputBEDChr,
          fromVcfBEDChr);

      logger.log(Level.INFO, "Adding padding to VCF regions");
      List<GenomeRange> paddedVcfBEDChr = generatePaddedGenomeRanges(filteredVcfBEDChr);

      logger.log(Level.INFO, "Merging overlapping VCF regions");
      List<GenomeRange> mergedVcfBEDChr = mergeOverlaps(paddedVcfBEDChr);

      logger.log(Level.INFO, String
          .format("Merging query regions with regions from VCF (%d records) from chromosome %s",
              mergedVcfBEDChr.size(), chromosome));
      List<GenomeRange> intermediateBED = mergeRegionsFromQueryBEDAndVariants(
          inputBEDChr, mergedVcfBEDChr, bedWindowSize);
      queryBED.addAll(massageBED(intermediateBED));
    }
    return queryBED;
  }
}
