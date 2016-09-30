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

import static com.google.common.base.Preconditions.checkArgument;

import com.google.cloud.genomics.utils.grpc.VariantUtils;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.genomics.v1.Range;
import com.google.genomics.v1.Variant;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange;
import com.verily.genomewarp.HomologousRangeOuterClass.TransformedVariants;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.RegionType;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.TargetStrand;
import com.verily.genomewarp.HomologousRangeOuterClass.TransformedVariants.TransformationStatus;
import com.verily.genomewarp.utils.Fasta;
import com.verily.genomewarp.utils.GenomeWarpUtils;
import htsjdk.samtools.util.SequenceUtil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.Set;

/**
 * Transforms variants in a query genome assembly to the corresponding representation in a target
 * genome assembly.
 *
 * <p>This class operates on variants within a single contiguous genome range on a query assembly.
 * Using the position and genotype information of those variants, along with the coordinates of the
 * genome region in both the query and target genome assemblies and the associated reference genome
 * sequences, it transforms the query variants into the representation of the same haplotypes on the
 * target genome.
 */
public final class TransformVariantsInRange {

  public static final Logger logger = Logger.getLogger(TransformVariantsInRange.class.getName());

  // Sentinel value for regions that cannot be successfully transformed.
  private static final TransformedVariants UNSUPPORTED_REGION =
      TransformedVariants.newBuilder().setStatus(TransformationStatus.FAILURE).build();

  // Sentinel value for a list of TransformationUnits that are not valid.
  private static final List<TransformationUnit> INVALID_TRANSFORMATION_UNITS = null;

  // Sentinel value for an allele that repeats beyond the boundary of the genome assembly region.
  private static final int ALLELE_REPEAT_END_NOT_FOUND = -1;

  /**
   * Takes a list of variants which all fall within the given region and performs
   * a variant transform, mapping them from one genome to another.
   *
   * @param region the region the input variants fall within
   * @param variants the list of variants to transform
   * @param callSetNames the callset names (sample names) of the input variants
   * @param queryFasta a fasta structure encoding the query genome
   * @param targetFasta a fasta structure encoding the target genome
   * @return a list of successfully transformed variants or null if the
   * variants cannot be transformed.
   */
  public static List<Variant> transformQueryToTargetVariants(HomologousRange region,
      List<Variant> variants, List<String> callSetNames, Fasta queryFasta, Fasta targetFasta) {
    List<Variant> successfulTransforms = null;

    TransformedVariants tVar = TransformVariantsInRange.transform(region, variants, callSetNames,
        queryFasta, targetFasta);

    if (tVar.getStatus() == TransformationStatus.SUCCESS) {
      // Successfully transfomed variants
      successfulTransforms = tVar.getVariantsList();
    }

    return successfulTransforms;
  }

  /**
   * Returns the variants present in the target genome assembly.
   *
   * <p>This function currently returns a {@code Pair} with a boolean value that indicates whether
   * the function operating on this type of range is valid. This is so that we can support the
   * staged rollout of different levels of transformation support--if the boolean is false then the
   * range and its (empty) variant list should not be considered a confidently-called region.
   *
   * @param homologousRange the pair of query and target assembly ranges to transform
   * @param queryVariants the list of Variants present in the query assembly
   * @param callSetName the name that should be given to each output {@code VariantCall}
   * @param queryFasta the {@code Fasta} for the query assembly
   * @param targetFasta the {@code Fasta} for the target assembly
   * @return a pair of Boolean (true if the region was successfully analyzed) and the list of
   *         {@code Variant}s present in the target assembly in this region.
   */
  public static TransformedVariants transform(HomologousRange homologousRange,
      List<Variant> queryVariants, String callSetName, Fasta queryFasta, Fasta targetFasta) {
    return transform(homologousRange, queryVariants, Arrays.asList(callSetName), queryFasta,
        targetFasta);
  }

  public static TransformedVariants transform(
      HomologousRange homologousRange,
      List<Variant> queryVariants,
      List<String> callSetNames,
      Fasta queryFasta,
      Fasta targetFasta) {
    validateInput(homologousRange, queryVariants);

    // Alignment-requiring regions are not yet supported.
    if (homologousRange.getRegionType() == RegionType.ALIGNMENT_REQUIRED) {
      return UNSUPPORTED_REGION;
    }

    // Multinucleotide variants are only supported if there are no reference genome changes and
    // they occur on the positive strand. Single nucleotide changes are supported in all cases.
    boolean hasOnlySnvs = hasOnlySingleNucleotideVariants(queryVariants);
    if (!hasOnlySnvs
        && (homologousRange.getRegionType() == RegionType.MISMATCHED_BASES
            || homologousRange.getTargetStrand() == TargetStrand.NEGATIVE_STRAND)) {
      return UNSUPPORTED_REGION;
    }

    List<ReferenceGenomeDifference> referenceGenomeChanges;
    // Optimization to avoid instantiating ReferenceGenomeDifference objects for regions and
    // variants in which we know for certain that reference genome changes will not occur.
    if (hasOnlySnvs && homologousRange.getRegionType() == RegionType.IDENTICAL) {
      referenceGenomeChanges = ImmutableList.of();
    } else {
      // Collect all reference genome changes induced by the coordinate liftover.
      List<ReferenceGenomeDifference> referenceGenomeChangesFromAssemblies =
          getReferenceChangesFromAssemblies(homologousRange, queryFasta, targetFasta);
      // Add all reference genome changes that were not properly left-shifted by the chain file.
      List<ReferenceGenomeDifference> referenceGenomeChangesFromVariants =
          getReferenceChangesFromVariants(
              homologousRange, queryVariants, queryFasta, targetFasta);
      referenceGenomeChanges =
          ImmutableList.copyOf(
              Iterables.concat(
                  referenceGenomeChangesFromAssemblies, referenceGenomeChangesFromVariants));
    }

    List<TransformationUnit> units =
        createTransformationUnits(referenceGenomeChanges, queryVariants, homologousRange);
    TransformedVariants targetVariants =
        getTargetVariantsForAllTransformationUnitsInRange(units, callSetNames);
    return returnValidOutput(homologousRange, targetVariants);
  }

  /**
   * Returns a list of all reference genome sequence differences present in the genome assemblies.
   *
   * @param homologousRange the pair of query and target assembly ranges of interest
   * @param queryFasta a reference to the query genome assembly
   * @param targetFasta a reference to the target genome assembly
   * @return a list of all reference genome sequence differences present in the genome assemblies
   */
  @VisibleForTesting
  static List<ReferenceGenomeDifference> getReferenceChangesFromAssemblies(
      HomologousRange homologousRange, Fasta queryFasta, Fasta targetFasta) {
    checkArgument(
        homologousRange.getRegionType() == RegionType.IDENTICAL
            || homologousRange.getRegionType() == RegionType.MISMATCHED_BASES,
        "Cannot identify reference changes from assembly regions of different sizes.");
    // Short circuit if this region has already been classified as having identical base pairs.
    if (homologousRange.getRegionType() == RegionType.IDENTICAL) {
      return ImmutableList.of();
    }

    ImmutableList.Builder<ReferenceGenomeDifference> differences = ImmutableList.builder();
    String queryDna =
        queryFasta.get(
            homologousRange.getQueryRange().getReferenceName(),
            homologousRange.getQueryRange().getStart(),
            homologousRange.getQueryRange().getEnd());
    String targetDna =
        targetFasta.get(
            homologousRange.getTargetRange().getReferenceName(),
            homologousRange.getTargetRange().getStart(),
            homologousRange.getTargetRange().getEnd());
    checkArgument(
        GenomeWarpUtils.isValidDna(queryDna),
        "Invalid query region DNA in %s: %s",
        homologousRange,
        queryDna);
    checkArgument(
        GenomeWarpUtils.isValidDna(targetDna),
        "Invalid target region DNA in %s: %s",
        homologousRange,
        targetDna);

    if (homologousRange.getTargetStrand() == TargetStrand.NEGATIVE_STRAND) {
      targetDna = SequenceUtil.reverseComplement(targetDna);
    }
    for (int i = 0; i < queryDna.length(); i++) {
      if (queryDna.charAt(i) != targetDna.charAt(i)) {
        differences.add(
            ReferenceGenomeDifference.create(
                homologousRange.getQueryRange().getStart() + i,
                queryDna.substring(i, i + 1),
                targetDna.substring(i, i + 1)));
      }
    }
    return differences.build();
  }

  private static List<ReferenceGenomeDifference> getReferenceChangesFromVariants(
      HomologousRange homologousRange,
      List<Variant> queryVariants,
      Fasta queryFasta,
      Fasta targetFasta) {
    // Short circuit if there are only SNVs, which cannot induce reference genome changes (see
    // http://go/genomewarp for details).
    if (hasOnlySingleNucleotideVariants(queryVariants)) {
      return ImmutableList.of();
    }

    checkArgument(
        homologousRange.getRegionType() == RegionType.IDENTICAL,
        "Finding reference changes from variants in different reference regions is unsupported.");
    checkArgument(
        homologousRange.getTargetStrand() == TargetStrand.POSITIVE_STRAND,
        "Finding reference changes from variants in reverse complement regions is unsupported.");

    ImmutableList.Builder<ReferenceGenomeDifference> differences = ImmutableList.builder();

    // Add the reference genome changes that occur outside the confidently-called region but, when
    // left-shifted to their VCF representation, cause a change within the confidently-called
    // region. These can only occur for indel variants.
    for (Variant queryVariant : queryVariants) {
      if (VariantUtils.IS_MULTI_NUCLEOTIDE.apply(queryVariant)) {
        ReferenceGenomeDifference leftShiftedChange =
            getReferenceChangeFromIndelInPositiveStrandIdenticalRegion(
                homologousRange, queryVariant, queryFasta, targetFasta);
        if (leftShiftedChange != ReferenceGenomeDifference.NO_DIFFERENCE) {
          differences.add(leftShiftedChange);
        }
      }
    }

    return differences.build();
  }

  /**
   * Returns a list of {@code TransformationUnit}s linking assembly changes to variant calls.
   *
   * <p>{@code TransformationUnit}s are the core object used to resolve changes between variants on
   * different genome assemblies. This function links reference genome assembly changes with
   * called variants in the same region.
   *
   * @param referenceGenomeDifferences reference genome assembly changes within a genomic region
   * @param variants a list of {@code Variant}s in the query assembly within the same genomic region
   * @return a list of {@code TransformationUnit}s that represent all variation in the region
   */
  @VisibleForTesting
  static List<TransformationUnit> createTransformationUnits(
      List<ReferenceGenomeDifference> referenceGenomeDifferences,
      List<Variant> variants,
      HomologousRange range) {
    ImmutableList.Builder<TransformationUnit> units = ImmutableList.builder();
    Set<Variant> unassignedVariants = new HashSet<>(variants);
    checkArgument(
        unassignedVariants.size() == variants.size(),
        "All variants in a single homologous range must be unique");
    for (ReferenceGenomeDifference referenceGenomeDifference : referenceGenomeDifferences) {
      ImmutableList.Builder<Variant> associatedVariants = ImmutableList.builder();
      for (Variant variant : variants) {
        if (referenceGenomeDifference.overlapsVariant(variant)) {
          associatedVariants.add(variant);
          // Try to remove this variant from the unassigned set. If it has already been removed,
          // it is attempting to be assigned to multiple reference genome changes. In this case we
          // don't know how to handle the region and give up.
          if (!unassignedVariants.remove(variant)) {
            logger.log(Level.WARNING,
                String.format("Complex region encountered in area with %s", variant));
            return INVALID_TRANSFORMATION_UNITS;
          }
        }
      }
      units.add(
          TransformationUnit.create(referenceGenomeDifference, associatedVariants.build(), range));
    }

    if (unassignedVariants.size() > 0) {
      units.add(
          TransformationUnit.create(
              ReferenceGenomeDifference.NO_DIFFERENCE,
              ImmutableList.copyOf(unassignedVariants),
              range));
    }

    return units.build();
  }

  /**
   * Returns all target genome assembly variants present in the list of {@code TransformationUnit}s.
   *
   * <p>The input {@code TransformationUnit} list must already have been checked for its internal
   * consistency (i.e. there are no overlapping units). If the list is invalid the region is
   * omitted.
   *
   * @param units the list of all {@code TransformationUnit}s in a single {@code HomologousRange}
   * @param callSetNames the callSetNames to use for the resulting variants
   * @return a pair of True (indicating that this region contains valid variants) and the list of
   *         {@code Variant}s present in the target assembly in this region.
   */
  @VisibleForTesting
  static TransformedVariants getTargetVariantsForAllTransformationUnitsInRange(
      List<TransformationUnit> units, List<String> callSetNames) {
    if (units == INVALID_TRANSFORMATION_UNITS) {
      return UNSUPPORTED_REGION;
    }

    List<Variant> allTargetVariants = new ArrayList<>();
    for (TransformationUnit unit : units) {
      List<Variant> targetVariants = unit.getTargetAssemblyVariants(callSetNames);
      if (targetVariants == TransformationUnit.INVALID_TARGET_VARIANTS) {
        // If any unit cannot be properly transformed, we must invalidate the entire region.
        return UNSUPPORTED_REGION;
      } else {
        allTargetVariants.addAll(targetVariants);
      }
    }

    // Sort output variants by chromosomal position to ensure consistent map shard outputs.
    Collections.sort(allTargetVariants, VariantUtils.CHROMOSOMAL_ORDER);

    return TransformedVariants.newBuilder()
        .setStatus(TransformationStatus.SUCCESS)
        .addAllVariants(allTargetVariants)
        .build();
  }

  /**
   * Returns the left-aligned reference genome difference within the region related to this indel.
   *
   * <p>A multinucleotide variant within a seemingly-identical genome region may have different copy
   * number in the query and target genome. This can occur because duplications of alleles can
   * extend past the boundary of a confident region.
   *
   * @param range the pair of query and target assembly ranges to transform
   * @param multiVariant the multi-nucleotide {@code Variant} present in the query assembly
   * @return the reference genome change occurring at the multiVariant position after left alignment
   */
  @VisibleForTesting
  static ReferenceGenomeDifference getReferenceChangeFromIndelInPositiveStrandIdenticalRegion(
      HomologousRange range,
      Variant multiVariant,
      Fasta queryFasta,
      Fasta targetFasta) {
    List<String> nonAnchorIndelAlleles = getNonAnchorIndelAlleles(multiVariant);

    // Complex variants can be trivially transformed, and will result in zero non-anchor indel
    // alleles being generated in this list.
    if (nonAnchorIndelAlleles.isEmpty()) {
      return ReferenceGenomeDifference.NO_DIFFERENCE;
    }

    long offsetFromRegionStart = multiVariant.getStart() - range.getQueryRange().getStart();
    int confidentRegionSize = (int) (range.getQueryRange().getEnd() - multiVariant.getStart());
    int querySize = confidentRegionSize;
    for (String allele : nonAnchorIndelAlleles) {
      if (allele.length() + 1 > querySize) {
        querySize = allele.length() + 2;
      }
    }
    String queryDNA =
        queryFasta.get(
            range.getQueryRange().getReferenceName(),
            multiVariant.getStart(),
            multiVariant.getStart() + querySize);

    // Find the first position in the query sequence where the DNA is not a copy of any of the
    // variant alleles.
    int queryMaxNonAlleleOffset = 0;
    String maxAllele = null;
    for (String allele : nonAnchorIndelAlleles) {
      int alleleNonAlleleOffset =
          firstNonAlleleOffset(
              allele, range.getQueryRange(), multiVariant.getStart(), queryFasta, queryDNA);
      if (alleleNonAlleleOffset > queryMaxNonAlleleOffset) {
        queryMaxNonAlleleOffset = alleleNonAlleleOffset;
        maxAllele = allele;
      }
    }

    // Fast case where we know the DNA is identical already.
    if (queryMaxNonAlleleOffset < confidentRegionSize) {
      return ReferenceGenomeDifference.NO_DIFFERENCE;
    }

    // Otherwise, find the offset in the target reference genome.
    String targetDNA =
        targetFasta.get(
            range.getTargetRange().getReferenceName(),
            range.getTargetRange().getStart() + offsetFromRegionStart,
            range.getTargetRange().getStart() + offsetFromRegionStart + queryMaxNonAlleleOffset);

    int targetMaxNonAlleleOffset =
        firstNonAlleleOffset(
            maxAllele,
            range.getTargetRange(),
            range.getTargetRange().getStart() + offsetFromRegionStart,
            targetFasta,
            targetDNA);

    int numDifferentAlleleCopies =
        (queryMaxNonAlleleOffset - targetMaxNonAlleleOffset) / maxAllele.length();
    if (numDifferentAlleleCopies == 0) {
      return ReferenceGenomeDifference.NO_DIFFERENCE;
    } else {
      String shortAllele = multiVariant.getReferenceBases().substring(0, 1);
      int numCopies = Math.abs(numDifferentAlleleCopies);
      StringBuilder longAllele = new StringBuilder(1 + numCopies * maxAllele.length());
      longAllele.append(shortAllele);
      for (int i = 0; i < numCopies; i++) {
        longAllele.append(maxAllele);
      }

      if (numDifferentAlleleCopies < 0) {
        return ReferenceGenomeDifference.create(
            multiVariant.getStart(), shortAllele, longAllele.toString());
      } else {
        return ReferenceGenomeDifference.create(
            multiVariant.getStart(), longAllele.toString(), shortAllele);
      }
    }
  }

  /**
   * Returns the offset from the initial position at which this allele is no longer replicated.
   *
   * <p>This is a wrapper function around firstNonAlleleOffsetHelper, which does the heavy lifting.
   * This wrapper is used to avoid many calls to the expensive Fasta.get function.
   *
   * @param allele the allele of interest
   * @param range the entire span of the region enclosing the variant in the reference of interest
   * @param variantStart the position of the variant of interest on the assembly
   * @param fasta a Fasta object for the genome assembly
   * @param dna a cached string holding the dna from the initial genome position of the variant to
   *      the end of the confidently-called region in the query assembly
   * @return the first index from the variant at which the reference genome no longer matches the
   *      allele of interest
   */
  private static int firstNonAlleleOffset(
      String allele, Range range, long variantStart, Fasta fasta, String dna) {
    int offset = firstNonAlleleOffsetHelper(dna, allele);
    long chromosomeSize = fasta.getChromosomeSize(range.getReferenceName());
    long end = 0;
    while (offset == ALLELE_REPEAT_END_NOT_FOUND) {
      if (end == chromosomeSize) {
        // The DNA repeats all the way to the end of the chromosome.
        // Use the length of the remaining chromosome as the offset.
        // This is a safe cast since the length of all chromosomes easily fit within an int.
        return (int) (end - variantStart);
      }
      end = Math.min(variantStart + (dna.length() * 2), chromosomeSize);
      String newQuery = fasta.get(range.getReferenceName(), variantStart, end);
      dna = newQuery;
      offset = firstNonAlleleOffsetHelper(dna, allele);
    }
    return offset;
  }

  /**
   * Returns the offset from the initial position at which this allele is no longer replicated, or
   * ALLELE_REPEAT_END_NOT_FOUND if the query DNA is a series of duplications of the allele.
   *
   * @param queryDNA the reference genome data to scan for duplicates of the allele
   * @param allele the allele to check duplications of
   * @return the first index in queryDNA for which it is no longer a replicate of allele
   */
  private static int firstNonAlleleOffsetHelper(String queryDNA, String allele) {
    // The counter starts at 1 because VCF conventions begin with an anchor base pair that is not
    // part of the indel variant.
    for (int queryIx = 1; queryIx < queryDNA.length(); queryIx++) {
      int alleleIx = (queryIx - 1) % allele.length();
      if (queryDNA.charAt(queryIx) != allele.charAt(alleleIx)) {
        return queryIx;
      }
    }
    return ALLELE_REPEAT_END_NOT_FOUND;
  }

  /**
   * Returns a list of all non-anchor indel alleles in a multi-nucleotide variant.
   *
   * <p>The VCF representation of a multi-nucleotide variant always includes an anchoring base.
   * For instance, an insertion of "T" at a reference base "C" is encoded as "C" --> "CT". This
   * function returns all of the non-anchored alleles of positive length (in this example, "T").
   *
   * <p>Complex alleles may not have the same anchor base as the reference. Those alleles can be
   * ignored since the anchor is guaranteed to be within the confidently-called region.
   *
   * @param variant the multi-nucleotide variant to examine
   * @return a list of all non-anchor alleles in a multi-nucleotide variant
   */
  private static List<String> getNonAnchorIndelAlleles(Variant variant) {
    ImmutableList.Builder<String> alleles = ImmutableList.builder();
    String reference = variant.getReferenceBases();
    if (reference.length() > 1) {
      alleles.add(reference.substring(1));
    }
    for (int i = 0; i < variant.getAlternateBasesCount(); i++) {
      String alt = variant.getAlternateBases(i);
      // Include indel alleles--those where the anchor base is the same as the reference.
      if (alt.length() > 1 && reference.charAt(0) == alt.charAt(0)) {
        alleles.add(alt.substring(1));
      }
    }
    return alleles.build();
  }

  /**
   * Ensures that the input range and variants are in a valid configuration.
   *
   * @param homologousRange the homologous range linking the query and target assemblies
   * @param queryVariants the list of variants within the homologous range on the query assembly.
   * @throws IllegalArgumentException if the homologous range is invalid or any query variant is not
   *         within the range.
   */
  private static void validateInput(HomologousRange homologousRange, List<Variant> queryVariants) {
    // Validate inputs.
    checkArgument(!homologousRange.getQueryRange().getReferenceName().isEmpty());
    checkArgument(!homologousRange.getTargetRange().getReferenceName().isEmpty());
    checkArgument(homologousRange.getRegionType() != RegionType.UNKNOWN_REGION_TYPE);
    checkArgument(homologousRange.getTargetStrand() != TargetStrand.UNKNOWN_TARGET_STRAND);

    for (Variant queryVariant : queryVariants) {
      // Sanity check all variants.
      checkArgument(
          (homologousRange
                  .getQueryRange()
                  .getReferenceName()
                  .equals(queryVariant.getReferenceName())
              && homologousRange.getQueryRange().getStart() <= queryVariant.getStart()
              && queryVariant.getStart() < homologousRange.getQueryRange().getEnd()),
          "Query variant is not within query range.");
    }
  }

  /**
   * Returns the valid target variants, or an empty list if the variants are not valid.
   *
   * <p>The transformations from query to target assembly generate a list of candidate target
   * assembly variants. Depending on the surrounding genome context, some indel variants on the
   * negative strand can have their position migrated, possibly outside the boundary of the
   * confidently-called region they are supposed to be enclosed within. This function checks that
   * all reported target variants lie within the appropriate confidently-called region, and if so
   * returns them. If any variant has migrated out of the region, the region is invalid and flagged
   * as such.
   *
   * @param homologousRange the homologous range linking the query and target assemblies
   * @param candidateOutput the candidate output to be analyzed
   * @return a pair of True (indicating that the output variants are valid) and the list of Variants
   *         present in the target assembly in this region.
   */
  private static TransformedVariants returnValidOutput(
      HomologousRange homologousRange, TransformedVariants candidateOutput) {
    if (candidateOutput.getStatus() != TransformationStatus.SUCCESS) {
      return candidateOutput;
    }

    for (Variant targetVariant : candidateOutput.getVariantsList()) {
      if (!(homologousRange
              .getTargetRange()
              .getReferenceName()
              .equals(targetVariant.getReferenceName())
          && homologousRange.getTargetRange().getStart() <= targetVariant.getStart()
          && targetVariant.getStart() < homologousRange.getTargetRange().getEnd())) {
        logger.log(Level.WARNING,
            String.format("Homologous range %s generates invalid variants", homologousRange));
        return UNSUPPORTED_REGION;
      }
    }
    return candidateOutput;
  }

  /** Returns True if and only if the list of variants is entirely composed of SNVs. */
  private static boolean hasOnlySingleNucleotideVariants(List<Variant> variants) {
    for (Variant variant : variants) {
      if (VariantUtils.IS_MULTI_NUCLEOTIDE.apply(variant)) {
        return false;
      }
    }
    return true;
  }

  private TransformVariantsInRange() {}
}
