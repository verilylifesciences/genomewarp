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
import com.google.genomics.v1.Variant;

/**
 * Container class for variation between the query and target reference genome assemblies.
 *
 * <p>Variation between reference assemblies influences the output of GenomeWarp by introducing or
 * eliminating variants present in the query assembly. This container class represents variation
 * between reference assemblies and the position of the variant in the query assembly. It omits
 * other fields that could be present (e.g. query chromosome, target chromosome, target position)
 * because these fields are inferred from other objects that are used in conjunction with this
 * class.
 */
@AutoValue
abstract class ReferenceGenomeDifference {
  // The zero-based inclusive position of the first base in the query assembly DNA sequence
  abstract long queryPosition();

  // The query assembly DNA sequence
  abstract String queryDna();

  // The target assembly DNA sequence
  abstract String targetDna();

  // Sentinel value for reference genome regions that are identical across assemblies.
  static final ReferenceGenomeDifference NO_DIFFERENCE =
      new AutoValue_ReferenceGenomeDifference(-1L, "X", "X");

  /**
   * Returns a {@code ReferenceGenomeDifference} object with the given parameters.
   *
   * @param position the zero-based inclusive position of the first base in the query sequence
   * @param query the query assembly DNA sequence
   * @param target the target assembly DNA sequence
   * @return a {@code ReferenceGenomeDifference} object with the given parameters
   */
  static ReferenceGenomeDifference create(long position, String query, String target) {
    // Validate the input fields. Indels are represented analogously to VCF format, where there is
    // an anchor base and then the variation.
    checkArgument(!query.isEmpty());
    checkArgument(!target.isEmpty());
    checkArgument(
        (query.length() == 1) || (target.length() == 1), "Complex changes are not permitted.");

    if (query.length() == target.length()) {
      // SNVs must have different base pairs.
      checkArgument(!query.toUpperCase().equals(target.toUpperCase()));
    } else {
      // Indels must have an identical anchor base.
      checkArgument(query.charAt(0) == target.charAt(0));
    }
    return new AutoValue_ReferenceGenomeDifference(position, query, target);
  }

  /** Returns true if and only if the reference change is a single nucleotide variant. */
  boolean isSnv() {
    return (!NO_DIFFERENCE.equals(this) && queryDna().length() == 1 && targetDna().length() == 1);
  }

  /** Returns true if and only if the reference change is a deletion in the target assembly. */
  boolean isDeletion() {
    return (!NO_DIFFERENCE.equals(this) && queryDna().length() > targetDna().length());
  }

  /** Returns true if and only if the reference change is an insertion in the target assembly. */
  boolean isInsertion() {
    return (!NO_DIFFERENCE.equals(this) && queryDna().length() < targetDna().length());
  }

  /**
   * Returns true if this query sequence overlaps the position of the input variant.
   *
   * <p>This function does not take into account the chromosome on which the variant lies, but only
   * compares the query sequence position and the variant position. It is the caller's job to ensure
   * this object and the input variant are on the same chromosome already.
   *
   * @param variant a {@code Variant} on the query assembly with which to check overlap
   * @return true if and only if this overlaps the variant on the query assembly
   */
  boolean overlapsVariant(Variant variant) {
    // Both queryPosition() and the variant are represented by zero-based half-open coordinates.
    return (queryPosition() < variant.getEnd()
        && queryPosition() + queryDna().length() > variant.getStart());
  }
}
