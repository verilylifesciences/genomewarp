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

package com.verily.genomewarp.utils;

import com.google.common.collect.ComparisonChain;
import com.google.genomics.v1.Position;
import com.google.genomics.v1.Range;
import com.google.genomics.v1.Variant;
import com.google.genomics.v1.VariantCall;
import com.verily.genomewarp.GenomeWarpSerial;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.RegionType;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.TargetStrand;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Utility functions for performing GenomeWarp
 */
public class GenomeWarpUtils {

  // The size in base pairs of each bucket in which to bin HomologousRanges and Variants.
  // Given ~4M variants in a whole human genome (~3B base pairs) implies that each bucket will have
  // BUCKET_SORT_SIZE / (3e9 / 4e6) = 133 variants in it.
  private static final long BUCKET_SORT_SIZE = 100_000L;
  public static final Pattern dnaSequence = Pattern.compile("[ACTGactg]+");
  public static final int VARIANT_CONTEXT_SIZE = 5;

  public static void fail(Logger logger, String message) {
    logger.log(Level.SEVERE, message);
    throw new RuntimeException(message);
  }

  private enum GenotypeCategory {
    REF_HOMOZYGOUS, // Diploid
    HETEROZYGOUS, // Diploid
    NON_REF_HOMOZYGOUS, // Diploid
    REF, // Haploid
    NON_REF // Haploid
  }

  private static final GenotypeCategory classifyGenotype(List<Integer> genotype) {
    boolean isHaploid = genotype.size() < 2;
    boolean hasRef = false;
    boolean allRef = true;
    for (Integer allele : genotype) {
      if (allele == 0) {
        hasRef = true;
      } else {
        allRef = false;
      }
    }

    if (allRef) {
      // Reference homozygous iff {0,0}, haploid reference if {0}
      return isHaploid ? GenotypeCategory.REF : GenotypeCategory.REF_HOMOZYGOUS;
    } else if (hasRef) {
      return GenotypeCategory.HETEROZYGOUS;
    }
    // We group all alternative alleles into non-reference for QC purposes.
    return isHaploid ? GenotypeCategory.NON_REF : GenotypeCategory.NON_REF_HOMOZYGOUS;
  }

  public static boolean hasVariation(Variant variant) {
    for (VariantCall call : variant.getCallsList()) {
      GenotypeCategory category = classifyGenotype(call.getGenotypeList());
      if (category != GenotypeCategory.REF_HOMOZYGOUS && category != GenotypeCategory.REF) {
        return true;
      }
    }
    return false;
  }

  /**
   * Returns true if the input string is a valid DNA string,
   * which we take to mean only containing the characters ACTGactg.
   *
   * @param inString the input string sequence to test
   * @return a boolean indicating whether the input is a DNA sequence
   */
  public static boolean isValidDna(String inString) {
    final byte[] in = StringUtil.stringToBytes(inString);
    for (int i = 0; i < in.length; i++) {
      if (!SequenceUtil.isValidBase(in[i])) {
        return false;
      }
    }
    return true;
  }

  /**
   * Defines an ordering on variants by comparing on the chromosome then
   * on the start position of the variant.
   */
  public static class VariantComparator implements Comparator<Variant> {
    @Override
    public int compare(Variant varA, Variant varB) {
      return ComparisonChain.start()
          .compare(varA.getReferenceName(), varB.getReferenceName())
          .compare(varA.getStart(), varB.getStart())
          .compare(varA.getEnd(), varB.getEnd())
          .result();
    }
  }

  /**
   * Defines an ordering on Range protos by comparing them on the reference
   * name, the start, then the end
   */
  public static class RangeComparator implements Comparator<Range> {
    @Override
    public int compare(Range a, Range b) {
      return ComparisonChain.start()
          .compare(a.getReferenceName(), b.getReferenceName())
          .compare(a.getStart(), b.getStart())
          .compare(a.getEnd(), b.getEnd())
          .result();
    }
  }

  /**
   * Defines a ordering on HomologousRange protos by comparing them on the
   * query chromosome, then target chromosome
   */
  public static class HomologousRangeComparator implements Comparator<HomologousRange> {
    @Override
    public int compare(HomologousRange a, HomologousRange b) {
      return ComparisonChain.start()
          .compare(a.getQueryRange(), b.getQueryRange(), new RangeComparator())
          .compare(a.getTargetRange(), b.getTargetRange(), new RangeComparator())
          .result();
    }
  }

  /**
   * Takes a list of variants and returns a new list of variants,
   * equal but without the info field present.
   *
   * @param variants the input list of variants
   * @return a list of variants without the info field
   */
  public static List<Variant> removeInfoField(List<Variant> variants) {

    List<Variant> toReturn = new ArrayList<>();

    for (Variant currVariant : variants) {
      // Pulls out just the fields we need from verbose variant record
      toReturn.add(Variant.newBuilder()
          .setReferenceName(currVariant.getReferenceName())
          .setStart(currVariant.getStart())
          .setEnd(currVariant.getEnd())
          .addAllNames(currVariant.getNamesList())
          .setReferenceBases(currVariant.getReferenceBases())
          .addAllAlternateBases(currVariant.getAlternateBasesList())
          .setQuality(currVariant.getQuality())
          .addAllFilter(currVariant.getFilterList())
          .addAllCalls(currVariant.getCallsList())
          .build());
    }

    return toReturn;
  }

  /**
   * Returns a position proto encoding the input reference name
   * and which chunk of size BUCKET_SORT_SIZE the start belongs in.
   *
   * @param refName the reference name of the returned position
   * @param start the start coordinate used to decide the chunk index
   * @return position
   */
  public static Position getPosition(String refName, long start) {
    return Position.newBuilder().setReferenceName(refName)
        .setPosition(start / BUCKET_SORT_SIZE).build();
  }

  private static Map<Position, List<HomologousRange>> associateRegionsWithPosition(
      List<HomologousRange> regions) {
    Map<Position, List<HomologousRange>> positionToRegions = new HashMap<>();
    for (HomologousRange currRegion : regions) {
      String reference = currRegion.getQueryRange().getReferenceName();
      long startBucket = currRegion.getQueryRange().getStart() / BUCKET_SORT_SIZE;
      long endBucket = currRegion.getQueryRange().getEnd() / BUCKET_SORT_SIZE;
      for (long i = startBucket; i <= endBucket; i++) {
        Position bucket = Position.newBuilder().setReferenceName(reference).setPosition(i).build();

        List<HomologousRange> regionList = null;
        if (positionToRegions.containsKey(bucket)) {
          regionList = positionToRegions.get(bucket);
        } else {
          regionList = new ArrayList<>();
        }

        regionList.add(currRegion);
        positionToRegions.put(bucket, regionList);
      }
    }

    return positionToRegions;
  }


  /**
   * Given a BufferedReader encoding a file of variants, asssociate the string
   * representation of the variants with a region. This association is done by
   * mapping each region to all variants which share the same reference name
   * (chromosome) and who's start coordiante is within the region.
   *
   * @param regions the list of regions to use as the map's keys
   * @param vcfFile a BufferedReader encoding the input variants
   * @return a map from regions to variant strings
   */
  public static Map<HomologousRange, List<String>> associateVariantsWithRange(
      List<HomologousRange> regions, BufferedReader vcfFile) throws IOException {
    Map<HomologousRange, List<String>> toReturn = new HashMap<>();

    // First associate each Region with one or more Positions
    Map<Position, List<HomologousRange>> positionToRegions =
        associateRegionsWithPosition(regions);

    // Add initial empty lists
    for (HomologousRange currRegion : regions) {
      toReturn.put(currRegion, new ArrayList<String>());
    }

    // We iterate through the variants. For each variant we get the corresponding
    // bucket, then check through the Regions in the bucket that suit our variant
    //
    // Furthermore, note that regions without any variants are still included, since
    // we previously added every region with an initially empty list of variants to toReturn.
    String line;
    while ((line = vcfFile.readLine()) != null) {
      if (line.length() == 0 || line.charAt(0) == '#') {
        // Header or blank line
        continue;
      }
      String[] vcfArray = line.split("\t");
      String refName = vcfArray[0];
      Long start = 0L;
      try {
        start = Long.parseLong(vcfArray[1]) - 1;
      } catch (NumberFormatException ex) {
        continue;
      }
      List<HomologousRange> correspondingRegions =
          positionToRegions.get(getPosition(refName, start));

      if (correspondingRegions == null) {
        // Variant occurs outside the confidently-called regions. These cannot be accurately
        // transformed, and so are ignored.
        continue;
      }

      for (HomologousRange currRegion : correspondingRegions) {
        if (start >= currRegion.getQueryRange().getStart()
            && start < currRegion.getQueryRange().getEnd()) {

          List<String> vcfList = toReturn.get(currRegion);
          vcfList.add(line);
          toReturn.put(currRegion, vcfList);
        }
      }
    }

    return toReturn;
  }

  /**
   * Duplication of associateVariantsWithRange. Required as we need to be able to deal
   * with variants formatted as strings as well as in true variant format.
   *
   * @param regions the list of regions to use as the map's keys
   * @param variants the input variants formatted as variant protos.
   * @return a map from regions to variant protos
   */
  public static Map<HomologousRange, List<Variant>> associateVariantsWithRange(
      List<HomologousRange> regions, List<Variant> variants) {
    Map<HomologousRange, List<Variant>> toReturn = new HashMap<>();

    Map<Position, List<HomologousRange>> positionToRegions = associateRegionsWithPosition(regions);

    // Add initial empty lists
    for (HomologousRange currRegion : regions) {
      toReturn.put(currRegion, new ArrayList<Variant>());
    }

    for (Variant currVariant : variants) {
      String refName = currVariant.getReferenceName();
      Long start = currVariant.getStart();
      List<HomologousRange> correspondingRegions =
          positionToRegions.get(getPosition(refName, start));

      if (correspondingRegions == null) {
        continue;
      }

      for (HomologousRange currRegion : correspondingRegions) {
        if (start >= currRegion.getQueryRange().getStart()
            && start < currRegion.getQueryRange().getEnd()) {

          List<Variant> variantList = toReturn.get(currRegion);
          if (variantList == null) {
            variantList = new ArrayList<>();
          }
          variantList.add(currVariant);
          toReturn.put(currRegion, variantList);
        }
      }
    }

    return toReturn;
  }

  public static String homologousToString(HomologousRange range) {
    Range queryRange = range.getQueryRange();
    Range targetRange = range.getTargetRange();
    return queryRange.getReferenceName() + "\t" + queryRange.getStart() + "\t" + queryRange.getEnd()
        + "\t" + targetRange.getReferenceName() + "\t" + targetRange.getStart() + "\t"
        + targetRange.getEnd() + "\t"
        + (range.getTargetStrand() == TargetStrand.POSITIVE_STRAND ? "+" : "-") + "\t"
        + range.getRegionType().toString();
  }

  /**
   * Parses a HomologousRange proto from an input string. This string must be tab separated with
   * each column defined as follows.
   *
   * <p> Col. 1: The query reference name (chromosome) <br>
   * Col. 2: The start position of the query range (long) <br>
   * Col. 3: The end position of the query range (long) <br>
   * Col. 4: The target reference name (chromosome) <br>
   * Col. 5: The start position of the target range (long) <br>
   * Col. 6: The end position of the target range (long) <br>
   * Col. 7: The strand (either + or -) <br>
   * Col. 8: The region type, encoded as one of IDENTICAL, ALIGNMENT_REQUIRED, or
   * MISMATCHED_BASES
   *
   * @param range the string encoding of the range
   * @return a HomologousRange parsed from the string
   */
  public static HomologousRange homologousFromString(String range) {
    String[] parts = range.split("\t");
    if (parts.length != 8) {
      throw new IllegalArgumentException("input string must have 9 sections separated by tabs");
    }

    Long qStart = 0L, qEnd = 0L, tStart = 0L, tEnd = 0L;
    try {
      qStart = Long.parseLong(parts[1]);
      qEnd = Long.parseLong(parts[2]);
      tStart = Long.parseLong(parts[4]);
      tEnd = Long.parseLong(parts[5]);
    } catch (NumberFormatException ex) {
      throw new IllegalArgumentException("Region.fromString: " + ex.getMessage());
    }
    boolean strand = parts[6].equals("+");

    return HomologousRange.newBuilder()
        .setQueryRange(Range.newBuilder().setReferenceName(parts[0]).setStart(qStart).setEnd(qEnd))
        .setTargetRange(Range.newBuilder().setReferenceName(parts[3]).setStart(tStart).setEnd(tEnd))
        .setTargetStrand(strand ? TargetStrand.POSITIVE_STRAND : TargetStrand.NEGATIVE_STRAND)
        .setRegionType(RegionType.valueOf(parts[7]))
        .build();
  }

  public static<T> T nextOrNull(Iterator<T> it) {
    if (it.hasNext()) {
      return it.next();
    }
    return null;
  }

}
