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

import com.google.auto.value.AutoValue;
import com.google.cloud.genomics.utils.grpc.VariantUtils;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableList;
import com.google.genomics.v1.Range;
import com.google.genomics.v1.Variant;
import com.google.genomics.v1.VariantCall;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.RegionType;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.TargetStrand;
import htsjdk.samtools.util.SequenceUtil;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.TreeSet;

/**
 * The fundamental unit of variant changes across genome assemblies.
 *
 * <p>Each {@code TransformationUnit} links one reference genome change to zero or more
 * {@code Variant}s within the {@code HomologousRange} on the query genome assembly that it spans.
 * This information is used to translate those sequence differences into the list of {@code Variant}
 * objects it corresponds to on the target genome assembly.
 */
@AutoValue
abstract class TransformationUnit {

  public static final Logger logger = Logger.getLogger(TransformationUnit.class.getName());

  // The encoding used by google.genomics.v1.VariantCall.genotype to indicate a no-call.
  private static final int NO_CALL = -1;

  // The encoding used by google.genomics.v1.Variant to indicate a passing genotype call.
  private static final String PASS_FILTER = "PASS";

  /** Sentinel value indicating that the list of target genome Variants is invalid. */
  static final List<Variant> INVALID_TARGET_VARIANTS = null;

  /** The single reference genome difference to consider. */
  abstract ReferenceGenomeDifference referenceGenomeDifference();

  /** The list of zero or more query Variants affected by this referenceGenomeDifference. */
  abstract List<Variant> variants();

  /** The homologous range in which the variants and reference genome difference lie. */
  abstract HomologousRange range();

  /**
   * Returns a {@code TransformationUnit} object with the given parameters.
   *
   * @param referenceGenomeDifference the single {@code ReferenceGenomeDifference} object
   * @param variants a list of {@code Variant}s that overlap the referenceGenomeDifference
   * @param range the {@code HomologousRange} in which the referenceGenomeDifference and variants
   *        occur
   */
  static TransformationUnit create(
      ReferenceGenomeDifference referenceGenomeDifference,
      List<Variant> variants,
      HomologousRange range) {
    // Sanity check all input values.
    Range queryRange = range.getQueryRange();
    long queryStart = queryRange.getStart();
    long queryEnd = queryRange.getEnd();
    checkArgument(
        (referenceGenomeDifference == ReferenceGenomeDifference.NO_DIFFERENCE)
            || ((queryStart <= referenceGenomeDifference.queryPosition())
                && (referenceGenomeDifference.queryPosition() < queryEnd)),
        "The ReferenceGenomeDifference must lie within the query range");
    for (Variant variant : variants) {
      checkArgument(
          queryRange.getReferenceName().equals(variant.getReferenceName())
              && queryStart <= variant.getStart()
              && variant.getStart() < queryEnd,
          "Every variant must be inside the query homologous range");
      checkArgument(
          referenceGenomeDifference == ReferenceGenomeDifference.NO_DIFFERENCE
              || referenceGenomeDifference.overlapsVariant(variant),
          "Either the reference must be unchanged or all variants must overlap it");
    }
    return new AutoValue_TransformationUnit(referenceGenomeDifference, variants, range);
  }

  /**
   * Returns the list of target assembly {@code Variant}s induced by this transformation unit.
   *
   * <p>At the moment only a subset of possible combinations of reference genome changes, linked
   * variants, and homologous range types can produce a valid list of target assembly variants. For
   * all other combinations, the sentinel value INVALID_TARGET_VARIANTS is returned, indicating that
   * the results of this attempted transformation should be ignored by downstream processing
   * pipelines.
   *
   * @param callSetNames the names of the callsets to set all variants to
   */
  List<Variant> getTargetAssemblyVariants(List<String> callSetNames) {
    // We can only transform regions with identical sizes, which are those with region type
    // IDENTICAL or MISMATCHED_BASES. Other types are not yet supported.
    if (range().getRegionType() != RegionType.IDENTICAL
        && range().getRegionType() != RegionType.MISMATCHED_BASES) {
      logger.log(Level.WARNING,
          String.format("Cannot transform regions with non-identical sizes: %s", this));
      return INVALID_TARGET_VARIANTS;
    }
    // Transformation of multinucleotide variants on the negative strand requires left-shifting
    // and an anchor base change, which is not yet supported.
    if (hasMultinucleotideVariation()
        && range().getTargetStrand() == TargetStrand.NEGATIVE_STRAND) {
      logger.log(Level.WARNING,
          String.format("Cannot transform multinucleotide variants on - strand: %s", this));
      return INVALID_TARGET_VARIANTS;
    }
    // Similar to the above, if the reference genome is an indel and it aligns on the negative
    // strand, we need to perform left-shifting and an anchor base change, which is not yet
    // supported.
    if ((referenceGenomeDifference().isInsertion() || referenceGenomeDifference().isDeletion())
        && range().getTargetStrand() == TargetStrand.NEGATIVE_STRAND) {
      logger.log(Level.WARNING, String.format("Cannot transform multinucleotide reference changes "
          + "on - strand: %s", this));
      return INVALID_TARGET_VARIANTS;
    }
    // Complex changes with multiple variants for a single reference change require complex logic
    // that is not yet supported.
    boolean identicalReferences =
        referenceGenomeDifference() == ReferenceGenomeDifference.NO_DIFFERENCE;
    if (!identicalReferences && variants().size() > 1) {
      logger.log(Level.WARNING, String.format("Cannot transform variant reference changes with "
          + "multiple variants: %s", this));
      return INVALID_TARGET_VARIANTS;
    }

    List<Variant> retval;
    if (identicalReferences) {
      retval = getTargetAssemblyVariantsInUnchangedGenome();
    } else if (variants().isEmpty()) {
      retval = getTargetAssemblyVariantsFromReferenceChangeOnly();
    } else {
      // The reference has changed, and a single variant overlaps the reference change.
      checkArgument(variants().size() == 1);
      Variant variant = variants().get(0);
      if (referenceGenomeDifference().isSnv() && !VariantUtils.IS_MULTI_NUCLEOTIDE.apply(variant)) {
        retval = getTargetAssemblyVariantsFromDualSNVs();
      } else if (variant.getAlternateBasesCount() == 1
          && referenceGenomeDifference().queryDna().equals(variant.getReferenceBases())
          && referenceGenomeDifference().targetDna().equals(variant.getAlternateBases(0))
          && range().getTargetStrand() == TargetStrand.POSITIVE_STRAND) {
        retval = getTargetAssemblyVariantsFromMatchingPositiveStrandIndels();
      } else {
        logger.log(Level.WARNING, String.format("Cannot transform variant reference changes with "
            + "this single variant: %s", this));
        return INVALID_TARGET_VARIANTS;
      }
    }
    return addCallSetNameToVariants(retval, callSetNames);
  }

  /**
   * Returns the list of target assembly {@code Variant}s induced by this transformation unit.
   *
   * <p>This function should only be called if it satisfies the following criteria:
   *   <ul>
   *     <li>The referenceGenomeDifference() == NO_DIFFERENCE
   *     <li>Either all variants are single nucleotide variants, or the change occurs on the
   *         positive strand.
   *   </ul>
   */
  private List<Variant> getTargetAssemblyVariantsInUnchangedGenome() {
    // If the reference genomes are unchanged, and we don't need to worry about left-shifting and
    // changing anchor bases for negative strand multinucleotide variants, the transformation is
    // trivial.
    ImmutableList.Builder<Variant> builder = ImmutableList.builder();
    for (Variant queryVariant : variants()) {
      builder.add(
          simpleCoordinateTransform(
              range(),
              queryVariant,
              queryVariant.getStart(),
              queryVariant.getReferenceBases(),
              queryVariant.getAlternateBasesList()));
    }
    return builder.build();
  }

  /**
   * Returns the list of target assembly {@code Variant}s induced by this transformation unit.
   *
   * <p>This function should only be called if it satisfies the following criteria:
   *   <ul>
   *     <li>The variants() list is empty.
   *     <li>Either the referenceGenomeDifference() is an SNV, or the change occurs on the positive
   *         strand.
   *   </ul>
   */
  private List<Variant> getTargetAssemblyVariantsFromReferenceChangeOnly() {
    // The individual is homozygous for the query reference sequence. Create a variant at the target
    // position for which the individual is homozygous for the alternate allele.
    Variant template =
        Variant.newBuilder()
            .addCalls(VariantCall.newBuilder().addGenotype(1).addGenotype(1))
            .addFilter(PASS_FILTER)
            .build();

    return ImmutableList.of(
        simpleCoordinateTransform(
            range(),
            template,
            referenceGenomeDifference().queryPosition(),
            referenceGenomeDifference().targetDna(),
            ImmutableList.of(referenceGenomeDifference().queryDna())));
  }

  /**
   * Returns the list of target assembly {@code Variant}s induced by this transformation unit.
   *
   * <p>This function should only be called if it satisfies the following criteria:
   *   <ul>
   *     <li>The variants() list contains a single SNV.
   *     <li>The referenceGenomeDifference() is an SNV.
   *   </ul>
   */
  private List<Variant> getTargetAssemblyVariantsFromDualSNVs() {
    Variant variant = variants().get(0);
    checkArgument(
        referenceGenomeDifference().queryDna().equals(variant.getReferenceBases()),
        "Different reference for query variant and reference in %s",
        range());

    // Create a map keyed by the genotype index with the value the string allele of that genotype in
    // the query genome assembly
    Map<Integer, String> queryGenotypeIndexToBase = new HashMap<>();
    queryGenotypeIndexToBase.put(0, variant.getReferenceBases());
    for (int i = 0; i < variant.getAlternateBasesCount(); i++) {
      queryGenotypeIndexToBase.put(i + 1, variant.getAlternateBases(i));
    }

    // Identify the target genome reference and all alternate alleles
    String targetReferenceBase = referenceGenomeDifference().targetDna();
    TreeSet<String> targetAltAlleleSet = new TreeSet<>();
    // Don't need to add referenceGenomeDifference().queryDna() since we verified it equals
    // variant.getReferenceBases()
    targetAltAlleleSet.add(variant.getReferenceBases());
    targetAltAlleleSet.addAll(variant.getAlternateBasesList());
    targetAltAlleleSet.remove(targetReferenceBase);
    List<String> targetAltAlleles = ImmutableList.copyOf(targetAltAlleleSet);

    // Create a map from the string allele in the target assembly to its genotype index
    Map<String, Integer> targetBaseToGenotypeIndex = new HashMap<>();
    targetBaseToGenotypeIndex.put(targetReferenceBase, 0);
    int i = 1;
    for (String altAllele : targetAltAlleles) {
      targetBaseToGenotypeIndex.put(altAllele, i++);
    }

    List<VariantCall> callsList = variant.getCallsList();
    Variant.Builder variantTemplate = Variant.newBuilder(variant)
        .clearReferenceBases()
        .clearAlternateBases()
        .clearCalls();
    for (VariantCall currCall : callsList) {
      // We only need to clear the genotype field. In particular, genotype likelihoods and phaseset
      // can remain the same since we repopulate genotypes in the same order.
      VariantCall.Builder callBuilder = VariantCall.newBuilder(currCall).clearGenotype();
      for (i = 0; i < currCall.getGenotypeCount(); i++) {
        int queryGenotypeIndex = currCall.getGenotype(i);
        if (queryGenotypeIndex == NO_CALL) {
          callBuilder.addGenotype(queryGenotypeIndex);
        } else {
          int targetGenotypeIndex =
              targetBaseToGenotypeIndex.get(queryGenotypeIndexToBase.get(queryGenotypeIndex));
          callBuilder.addGenotype(targetGenotypeIndex);
        }
      }
      variantTemplate.addCalls(callBuilder);
    }

    return ImmutableList.of(
        simpleCoordinateTransform(
            range(), variantTemplate.build(), variant.getStart(), targetReferenceBase,
            targetAltAlleles));
  }

  /**
   * Returns the list of target assembly {@code Variant}s induced by this transformation unit.
   *
   * <p>This function should only be called if it satisfies the following criteria:
   *   <ul>
   *     <li>The variants() list contains a single indel.
   *     <li>The referenceGenomeDifference() is an indel.
   *     <li>The alleles of the variant and the reference genome difference are identical.
   *     <li>The target assembly maps to the positive strand.
   *   </ul>
   */
  private List<Variant> getTargetAssemblyVariantsFromMatchingPositiveStrandIndels() {
    Variant variant = variants().get(0);
    String targetReferenceBases = referenceGenomeDifference().targetDna();
    String targetAlternateBases = referenceGenomeDifference().queryDna();

    Variant.Builder template = Variant.newBuilder(variant).clearCalls();
    List<VariantCall> callsList = variant.getCallsList();
    for (VariantCall currCall : callsList) {
      // Genotype likelihoods and phaseset can remain the same
      // since we repopulate in the same order.
      VariantCall.Builder callBuilder = VariantCall.newBuilder(currCall).clearGenotype();
      for (int i = 0; i < currCall.getGenotypeCount(); i++) {
        int queryGenotypeIndex = currCall.getGenotype(i);
        if (queryGenotypeIndex == NO_CALL) {
          callBuilder.addGenotype(queryGenotypeIndex);
        } else {
          checkArgument(0 <= queryGenotypeIndex && queryGenotypeIndex <= 1);
          callBuilder.addGenotype(1 - queryGenotypeIndex);
        }
      }
      template.addCalls(callBuilder);
    }

    return ImmutableList.of(
        simpleCoordinateTransform(
            range(),
            template.build(),
            variant.getStart(),
            targetReferenceBases,
            ImmutableList.of(targetAlternateBases)));
  }

  /**
   * Returns the target assembly variant representing this query variant.
   *
   * <p>This function performs the coordinate transformation and sequence complementing necessary
   * for variants on either strand. However, it does not perform any left-shifting or changing of
   * anchor bases for indels on the negative strand.
   *
   * @param range the homologous range linking the query and target assemblies
   * @param template the variant template to use in creation of the target variant
   * @param queryStartPosition the start position of the variant on the query assembly
   * @param targetPositiveStrandReferenceAllele the reference allele in the target variant if
   *        it is on the positive strand
   * @param targetPositiveStrandAlternateAlleles the list of alternate alleles in the target variant
   *        if it is on the positive strand
   * @return the variant of interest on the target assembly
   */
  @VisibleForTesting
  static Variant simpleCoordinateTransform(
      HomologousRange range,
      Variant template,
      long queryStartPosition,
      String targetPositiveStrandReferenceAllele,
      List<String> targetPositiveStrandAlternateAlleles) {
    long convertedQueryStartPosition = positionConversion(range, queryStartPosition);
    switch (range.getTargetStrand()) {
      case POSITIVE_STRAND:
        return Variant.newBuilder(template)
            .setReferenceName(range.getTargetRange().getReferenceName())
            .setStart(convertedQueryStartPosition)
            .setEnd(convertedQueryStartPosition + targetPositiveStrandReferenceAllele.length())
            .setReferenceBases(targetPositiveStrandReferenceAllele)
            .clearAlternateBases()
            .addAllAlternateBases(targetPositiveStrandAlternateAlleles)
            .build();
      case NEGATIVE_STRAND:
        ImmutableList.Builder<String> revcompAlternateAlleles = ImmutableList.builder();
        for (String allele : targetPositiveStrandAlternateAlleles) {
          revcompAlternateAlleles.add(SequenceUtil.reverseComplement(allele));
        }
        return Variant.newBuilder(template)
            .setReferenceName(range.getTargetRange().getReferenceName())
            // Note that start and end are swapped on the negative strand.
            .setStart(convertedQueryStartPosition - targetPositiveStrandReferenceAllele.length())
            .setEnd(convertedQueryStartPosition)
            .setReferenceBases(SequenceUtil.reverseComplement(targetPositiveStrandReferenceAllele))
            .clearAlternateBases()
            .addAllAlternateBases(revcompAlternateAlleles.build())
            .build();
      case UNKNOWN_TARGET_STRAND:
        // Fall through to default exception.
      default:
        throw new IllegalArgumentException("Cannot transform without known strand: " + range);
    }
  }

  /** Converts a position in the query to its analogous position in the target. */
  private static long positionConversion(HomologousRange range, long queryPosition) {
    long queryStart = range.getQueryRange().getStart();
    long queryEnd = range.getQueryRange().getEnd();
    checkArgument(
        (queryStart <= queryPosition) && (queryEnd > queryPosition),
        "Invalid query position to convert: %s in [%s, %s)",
        queryPosition,
        queryStart,
        queryEnd);
    switch (range.getTargetStrand()) {
      case POSITIVE_STRAND:
        return range.getTargetRange().getStart() + (queryPosition - queryStart);
      case NEGATIVE_STRAND:
        return range.getTargetRange().getEnd() - (queryPosition - queryStart);
      case UNKNOWN_TARGET_STRAND:
        // Fall through to default exception.
      default:
        throw new IllegalArgumentException("Cannot convert without known strand: " + range);
    }
  }

  /** Returns True iff any variant in this is a multinucleotide variant. */
  private boolean hasMultinucleotideVariation() {
    for (Variant variant : variants()) {
      if (VariantUtils.IS_MULTI_NUCLEOTIDE.apply(variant)) {
        return true;
      }
    }
    return false;
  }

  /** Returns a list of {@code Variant}s with the callSetNames field populated. */
  private static List<Variant> addCallSetNameToVariants(
      List<Variant> variants, List<String> callSetNames) {
    ImmutableList.Builder<Variant> builder = ImmutableList.builder();
    for (Variant variant : variants) {
      Variant.Builder namedVariant = Variant.newBuilder(variant).clearCalls();
      List<VariantCall> callsList = variant.getCallsList();
      int i = 0;
      for (VariantCall currCall : callsList) {
        VariantCall namedCall =
            VariantCall.newBuilder(currCall).setCallSetName(callSetNames.get(i++)).build();
        namedVariant.addCalls(namedCall);
      }
      builder.add(namedVariant.build());
    }
    return builder.build();
  }
}
