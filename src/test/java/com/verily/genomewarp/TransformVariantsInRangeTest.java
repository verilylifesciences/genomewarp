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

import static com.verily.genomewarp.utils.GenomeWarpTestUtils.createSomeCalls;
import static com.verily.genomewarp.utils.GenomeWarpTestUtils.makeVariant;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.genomics.v1.Range;
import com.google.genomics.v1.Variant;
import com.verily.genomewarp.ReferenceGenomeDifference;
import com.verily.genomewarp.TransformVariantsInRange;
import com.verily.genomewarp.TransformationUnit;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange;
import com.verily.genomewarp.HomologousRangeOuterClass.TransformedVariants;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.RegionType;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.TargetStrand;
import com.verily.genomewarp.HomologousRangeOuterClass.TransformedVariants.TransformationStatus;
import com.verily.genomewarp.utils.Fasta;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.junit.Rule;
import org.junit.Test;
import org.junit.experimental.runners.Enclosed;
import org.junit.rules.ExpectedException;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameter;
import org.junit.runners.Parameterized.Parameters;

@RunWith(Enclosed.class)
public final class TransformVariantsInRangeTest {

  private static final String QUERY_REFERENCE_PATH = "query.fasta";

  private static final String TARGET_REFERENCE_PATH = "target.fasta";

  private static final List<String> CALLSET_NAME = Arrays.asList("MyCallsetName1",
      "MyCallsetName2", "MyCallsetName3");

  private static final HomologousRange LARGEST_CHR1_HOMOLOG =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr1").setStart(1).setEnd(40))
          .setTargetRange(Range.newBuilder().setReferenceName("chr1_same").setStart(11).setEnd(50))
          .setTargetStrand(TargetStrand.POSITIVE_STRAND)
          .setRegionType(RegionType.IDENTICAL)
          .build();

  private static final HomologousRange TRUNCATED_CHR1_HOMOLOG =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr1").setStart(1).setEnd(13))
          .setTargetRange(Range.newBuilder().setReferenceName("chr1_same").setStart(11).setEnd(23))
          .setTargetStrand(TargetStrand.POSITIVE_STRAND)
          .setRegionType(RegionType.IDENTICAL)
          .build();

  private static final HomologousRange LARGEST_CHR2_HOMOLOG =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr2").setStart(1).setEnd(74))
          .setTargetRange(Range.newBuilder().setReferenceName("chr2_same").setStart(11).setEnd(84))
          .setTargetStrand(TargetStrand.POSITIVE_STRAND)
          .setRegionType(RegionType.IDENTICAL)
          .build();

  private static final HomologousRange LARGEST_CHR2_REVCOMP_HOMOLOG =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr2").setStart(1).setEnd(43))
          .setTargetRange(
              Range.newBuilder().setReferenceName("chr2_revcomp").setStart(10).setEnd(52))
          .setTargetStrand(TargetStrand.NEGATIVE_STRAND)
          .setRegionType(RegionType.IDENTICAL)
          .build();

  private static final HomologousRange CHR2_MISMATCHED_BASES_HOMOLOG =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr2").setStart(0).setEnd(75))
          .setTargetRange(
              Range.newBuilder().setReferenceName("chr2_mismatched_bases").setStart(0).setEnd(75))
          .setTargetStrand(TargetStrand.POSITIVE_STRAND)
          .setRegionType(RegionType.MISMATCHED_BASES)
          .build();

  private static final HomologousRange CHR2_REVCOMP_MISMATCHED_BASES_HOMOLOG =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr2").setStart(0).setEnd(74))
          .setTargetRange(
              Range.newBuilder().setReferenceName("chr2_revcomp_mismatch").setStart(4).setEnd(78))
          .setTargetStrand(TargetStrand.NEGATIVE_STRAND)
          .setRegionType(RegionType.MISMATCHED_BASES)
          .build();

  private static final HomologousRange LARGEST_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr2").setStart(1).setEnd(19))
          .setTargetRange(
              Range.newBuilder().setReferenceName("chr2_CTG_deletion").setStart(11).setEnd(29))
          .setTargetStrand(TargetStrand.POSITIVE_STRAND)
          .setRegionType(RegionType.IDENTICAL)
          .build();

  private static final HomologousRange SMALL_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr2").setStart(1).setEnd(9))
          .setTargetRange(
              Range.newBuilder().setReferenceName("chr2_CTG_deletion").setStart(11).setEnd(19))
          .setTargetStrand(TargetStrand.POSITIVE_STRAND)
          .setRegionType(RegionType.IDENTICAL)
          .build();

  private static final HomologousRange CHR2_PART_DELETION_IDENTICAL_HOMOLOG =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr2").setStart(1).setEnd(21))
          .setTargetRange(
              Range.newBuilder().setReferenceName("chr2_part_deletion").setStart(11).setEnd(31))
          .setTargetStrand(TargetStrand.POSITIVE_STRAND)
          .setRegionType(RegionType.IDENTICAL)
          .build();

  private static final HomologousRange CHR2_CTG_DELETION_ALIGNMENT_HOMOLOG =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr2").setStart(1).setEnd(74))
          .setTargetRange(
              Range.newBuilder().setReferenceName("chr2_CTG_deletion").setStart(11).setEnd(81))
          .setTargetStrand(TargetStrand.POSITIVE_STRAND)
          .setRegionType(RegionType.ALIGNMENT_REQUIRED)
          .build();

  private static final HomologousRange LARGEST_CHR2_CTG_INSERTION_IDENTICAL_HOMOLOG =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr2").setStart(1).setEnd(22))
          .setTargetRange(
              Range.newBuilder().setReferenceName("chr2_CTG_insertion").setStart(11).setEnd(33))
          .setTargetStrand(TargetStrand.POSITIVE_STRAND)
          .setRegionType(RegionType.IDENTICAL)
          .build();

  private static final HomologousRange CHR2_CTG_INSERTION_ALIGNMENT_HOMOLOG =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr2").setStart(1).setEnd(74))
          .setTargetRange(
              Range.newBuilder().setReferenceName("chr2_CTG_insertion").setStart(11).setEnd(87))
          .setTargetStrand(TargetStrand.POSITIVE_STRAND)
          .setRegionType(RegionType.ALIGNMENT_REQUIRED)
          .build();

  private static final Variant SINGLE_NUCLEOTIDE_HET1 =
      makeVariant("chr1", 17, "A", "G", createSomeCalls(CALLSET_NAME, 0, 1, 1, 1, 1, 0));

  private static final Variant SINGLE_NUCLEOTIDE_HET2 =
      makeVariant("chr1", 21, "G", "T", createSomeCalls(CALLSET_NAME, 1, 0, 1, 0));

  private static final Variant INSERTION_HET1 =
      makeVariant("chr2", 1, "A", "ACTG", createSomeCalls(CALLSET_NAME, 0, 1));

  @RunWith(Parameterized.class)
  public static final class SuccessfulTransformTests {
    @Parameters(name = "{index}: successfulTransform({0}, {1})={2}")
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            // Various different identical transformations.
            {
              LARGEST_CHR1_HOMOLOG,
              Arrays.asList(
                  makeVariant("chr1", 3, "G", "T", createSomeCalls(CALLSET_NAME, 0, 1, 1, 0)),
                  makeVariant("chr1", 8, "C", "T", createSomeCalls(CALLSET_NAME, 2, 1, 2, 1))
                      .toBuilder()
                      .addAlternateBases("G")
                      .build()),
              Arrays.asList(
                  makeVariant("chr1_same", 13, "G", "T", createSomeCalls(CALLSET_NAME, 0, 1, 1, 0)),
                  makeVariant("chr1_same", 18, "C", "T", createSomeCalls(CALLSET_NAME, 2, 1, 2, 1))
                      .toBuilder()
                      .addAlternateBases("G")
                      .build())
            },
            {
              LARGEST_CHR1_HOMOLOG,
              Arrays.asList(
                  makeVariant("chr1", 8, "CG", "TG", createSomeCalls(CALLSET_NAME, 0, 2))
                      .toBuilder()
                      .addAlternateBases("C")
                      .build()),
              Arrays.asList(
                  makeVariant("chr1_same", 18, "CG", "TG", createSomeCalls(CALLSET_NAME, 0, 2))
                      .toBuilder()
                      .addAlternateBases("C")
                      .build())
            },
            {
              LARGEST_CHR1_HOMOLOG,
              Arrays.asList(
                  makeVariant("chr1", 8, "CG", "TG",
                    createSomeCalls(CALLSET_NAME, 0, 1, 0, 1, 0, 1))
                      .toBuilder()
                      .addAlternateBases("CGG")
                      .build()),
              Arrays.asList(
                  makeVariant("chr1_same", 18, "CG", "TG",
                    createSomeCalls(CALLSET_NAME, 0, 1, 0, 1, 0, 1))
                      .toBuilder()
                      .addAlternateBases("CGG")
                      .build())
            },
            {
              LARGEST_CHR1_HOMOLOG,
              Arrays.asList(makeVariant("chr1", 8, "C", "CT", createSomeCalls(CALLSET_NAME, 0, 1))),
              Arrays.asList(makeVariant("chr1_same", 18, "C", "CT",
                  createSomeCalls(CALLSET_NAME, 0, 1)))
            },
            {
              TRUNCATED_CHR1_HOMOLOG,
              Arrays.asList(
                  makeVariant("chr1", 8, "CG", "C", createSomeCalls(CALLSET_NAME, 1, 1)),
                  makeVariant("chr1", 10, "CA", "C", createSomeCalls(CALLSET_NAME, 0, 1))),
              Arrays.asList(
                  makeVariant("chr1_same", 18, "CG", "C", createSomeCalls(CALLSET_NAME, 1, 1)),
                  makeVariant("chr1_same", 20, "CA", "C", createSomeCalls(CALLSET_NAME, 0, 1)))
            },
            {
              LARGEST_CHR2_HOMOLOG,
              Arrays.asList(makeVariant("chr2", 1, "ACTG", "A",
                  createSomeCalls(CALLSET_NAME, -1, -1))),
              Arrays.asList(makeVariant("chr2_same", 11, "ACTG", "A",
                  createSomeCalls(CALLSET_NAME, -1, -1)))
            },
            {
              LARGEST_CHR2_HOMOLOG,
              Arrays.asList(makeVariant("chr2", 1, "A", "ACT",
                  createSomeCalls(CALLSET_NAME, 1, 1))),
              Arrays.asList(makeVariant("chr2_same", 11, "A", "ACT",
                  createSomeCalls(CALLSET_NAME, 1, 1)))
            },
            {
              LARGEST_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              Arrays.asList(makeVariant("chr2", 1, "ACT", "A",
                  createSomeCalls(CALLSET_NAME, 1, 1))),
              Arrays.asList(makeVariant("chr2_CTG_deletion", 11, "ACT", "A",
                  createSomeCalls(CALLSET_NAME, 1, 1)))
            },
            {
              LARGEST_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              Arrays.asList(makeVariant("chr2", 1, "A", "ACT",
                  createSomeCalls(CALLSET_NAME, 1, 0))),
              Arrays.asList(makeVariant("chr2_CTG_deletion", 11, "A", "ACT",
                  createSomeCalls(CALLSET_NAME, 1, 0)))
            },
            {
              LARGEST_CHR2_CTG_INSERTION_IDENTICAL_HOMOLOG,
              Arrays.asList(makeVariant("chr2", 1, "ACT", "A",
                  createSomeCalls(CALLSET_NAME, 0, 1, 1, 1))),
              Arrays.asList(
                  makeVariant("chr2_CTG_insertion", 11, "ACT", "A",
                  createSomeCalls(CALLSET_NAME, 0, 1, 1, 1)))
            },
            {
              CHR2_PART_DELETION_IDENTICAL_HOMOLOG,
              Arrays.asList(makeVariant("chr2", 1, "A", "ACTG",
                  createSomeCalls(CALLSET_NAME, 0, 1))),
              Arrays.asList(
                  makeVariant("chr2_part_deletion", 11, "A", "ACTG",
                  createSomeCalls(CALLSET_NAME, 0, 1)))
            },
            // Indels with edge effects from the ending of the confidently-called region.
            {
              LARGEST_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              Arrays.asList(
                  makeVariant("chr2", 1, "ACTG", "A", createSomeCalls(CALLSET_NAME, 1, 1)),
                  makeVariant("chr2", 10, "G", "A", createSomeCalls(CALLSET_NAME, 0, 1))),
              Arrays.asList(
                  makeVariant("chr2_CTG_deletion", 11, "A", "ACTG",
                    createSomeCalls(CALLSET_NAME, 0, 0)),
                  makeVariant("chr2_CTG_deletion", 20, "G", "A",
                    createSomeCalls(CALLSET_NAME, 0, 1)))
            }
          });
    }

    @Parameter(0)
    public HomologousRange testRange;

    @Parameter(1)
    public List<Variant> queryVariants;

    @Parameter(2)
    public List<Variant> expectedVariants;

    @Test
    public void supportedTransform() throws IOException {
      Fasta queryFasta = new Fasta(TransformVariantsInRangeTest.class.getClassLoader()
          .getResource(QUERY_REFERENCE_PATH).getFile());
      Fasta targetFasta = new Fasta(TransformVariantsInRangeTest.class.getClassLoader()
          .getResource(TARGET_REFERENCE_PATH).getFile());

      TransformedVariants result =
          TransformVariantsInRange.transform(
              testRange, queryVariants, CALLSET_NAME, queryFasta, targetFasta);
      assertTrue(result.getStatus().equals(TransformationStatus.SUCCESS));
      List<Variant> actualVariants = result.getVariantsList();
      assertEquals(actualVariants, expectedVariants);
      assertTrue(actualVariants.containsAll(expectedVariants)
          && expectedVariants.containsAll(actualVariants));
    }
  }

  @RunWith(Parameterized.class)
  public static final class UnsupportedTransformTests {
    @Parameters(name = "{index}: unsupportedTransform({0}, {1})")
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            // Alignment region.
            {
              CHR2_CTG_INSERTION_ALIGNMENT_HOMOLOG,
              Arrays.asList(makeVariant("chr2", 1, "A", "ACTG",
                  createSomeCalls(CALLSET_NAME, 0, 1)))
            },
            // Insertion on the negative strand.
            {
              LARGEST_CHR2_REVCOMP_HOMOLOG,
              Arrays.asList(makeVariant("chr2", 1, "A", "ACTG",
                  createSomeCalls(CALLSET_NAME, 0, 1)))
            },
            // Deletion with mismatched base pairs.
            {
              CHR2_MISMATCHED_BASES_HOMOLOG,
              Arrays.asList(makeVariant("chr2", 1, "AC", "A", createSomeCalls(CALLSET_NAME, 0, 1)))
            },
            // Indel causes multiple variants to map to a single genomic region change.
            {
              SMALL_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              Arrays.asList(
                  makeVariant("chr2", 1, "A", "ACTG", createSomeCalls(CALLSET_NAME, 0, 1)),
                  makeVariant("chr2", 2, "C", "T", createSomeCalls(CALLSET_NAME, 0, 1)))
            },
          });
    }

    @Parameter(0)
    public HomologousRange testRange;

    @Parameter(1)
    public List<Variant> queryVariants;

    @Test
    public void unsupportedTransform() throws IOException {
      Fasta queryFasta = new Fasta(TransformVariantsInRangeTest.class.getClassLoader()
          .getResource(QUERY_REFERENCE_PATH).getFile());
      Fasta targetFasta = new Fasta(TransformVariantsInRangeTest.class.getClassLoader()
          .getResource(TARGET_REFERENCE_PATH).getFile());

      TransformedVariants result =
          TransformVariantsInRange.transform(
              testRange, queryVariants, CALLSET_NAME, queryFasta, targetFasta);
      assertTrue(result.getStatus().equals(TransformationStatus.FAILURE));
      assertTrue(result.getVariantsList().size() == 0);
    }
  }

  @RunWith(Parameterized.class)
  public static final class ReferenceChangesFromAssembliesTests {
    @Parameters(name = "{index}: getReferenceChangesFromAssemblies({0})={1}")
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            {LARGEST_CHR1_HOMOLOG, Arrays.asList()},
            {LARGEST_CHR2_REVCOMP_HOMOLOG, Arrays.asList()},
            {
              CHR2_MISMATCHED_BASES_HOMOLOG,
              Arrays.asList(
                  ReferenceGenomeDifference.create(6, "T", "A"),
                  ReferenceGenomeDifference.create(31, "A", "T"),
                  ReferenceGenomeDifference.create(41, "A", "T"),
                  ReferenceGenomeDifference.create(54, "A", "T"),
                  ReferenceGenomeDifference.create(64, "A", "T"),
                  ReferenceGenomeDifference.create(73, "A", "T"))
            },
            {
              CHR2_REVCOMP_MISMATCHED_BASES_HOMOLOG,
              Arrays.asList(
                  ReferenceGenomeDifference.create(6, "T", "A"),
                  ReferenceGenomeDifference.create(31, "A", "T"),
                  ReferenceGenomeDifference.create(41, "A", "T"),
                  ReferenceGenomeDifference.create(54, "A", "T"),
                  ReferenceGenomeDifference.create(64, "A", "T"),
                  ReferenceGenomeDifference.create(73, "A", "T"))
            },
          });
    }

    @Parameter(0)
    public HomologousRange range;

    @Parameter(1)
    public List<ReferenceGenomeDifference> expectedRgds;

    @Test
    public void getReferenceChangesFromAssemblies() throws IOException {
      Fasta queryFasta = new Fasta(TransformVariantsInRangeTest.class.getClassLoader()
          .getResource(QUERY_REFERENCE_PATH).getFile());
      Fasta targetFasta = new Fasta(TransformVariantsInRangeTest.class.getClassLoader()
          .getResource(TARGET_REFERENCE_PATH).getFile());
      List<ReferenceGenomeDifference> actualRgds =
          TransformVariantsInRange.getReferenceChangesFromAssemblies(
              range, queryFasta, targetFasta);
      assertTrue(actualRgds.equals(expectedRgds));
    }
  }

  public static final class UnsupportedReferenceChangesFromAssembliesTests {

    @Rule public final ExpectedException thrown = ExpectedException.none();

    @Test
    public void unsupportedGetReferenceChangesFromAssemblies() throws IOException {
      thrown.expect(IllegalArgumentException.class);
      Fasta queryFasta = new Fasta(TransformVariantsInRangeTest.class.getClassLoader()
          .getResource(QUERY_REFERENCE_PATH).getFile());
      Fasta targetFasta = new Fasta(TransformVariantsInRangeTest.class.getClassLoader()
          .getResource(TARGET_REFERENCE_PATH).getFile());
      TransformVariantsInRange.getReferenceChangesFromAssemblies(
          CHR2_CTG_DELETION_ALIGNMENT_HOMOLOG, queryFasta, targetFasta);
    }
  }

  @RunWith(Parameterized.class)
  public static final class CreateTransformationUnitsTests {
    @Parameters(name = "{index}: createTransformationUnits({0}, {1})={2}")
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            {Arrays.asList(), Arrays.asList(), Arrays.asList()},
            // No reference genome changes.
            {
              Arrays.asList(),
              Arrays.asList(SINGLE_NUCLEOTIDE_HET1),
              Arrays.asList(
                  TransformationUnit.create(
                      ReferenceGenomeDifference.NO_DIFFERENCE,
                      Arrays.asList(SINGLE_NUCLEOTIDE_HET1),
                      LARGEST_CHR1_HOMOLOG))
            },
            {
              Arrays.asList(),
              Arrays.asList(SINGLE_NUCLEOTIDE_HET1, SINGLE_NUCLEOTIDE_HET2),
              Arrays.asList(
                  TransformationUnit.create(
                      ReferenceGenomeDifference.NO_DIFFERENCE,
                      Arrays.asList(SINGLE_NUCLEOTIDE_HET1, SINGLE_NUCLEOTIDE_HET2),
                      LARGEST_CHR1_HOMOLOG))
            },
            // Only reference genome changes.
            {
              Arrays.asList(ReferenceGenomeDifference.create(17, "A", "G")),
              Arrays.asList(),
              Arrays.asList(
                  TransformationUnit.create(
                      ReferenceGenomeDifference.create(17, "A", "G"),
                      Arrays.<Variant>asList(),
                      LARGEST_CHR1_HOMOLOG))
            },
            {
              Arrays.asList(
                  ReferenceGenomeDifference.create(17, "A", "G"),
                  ReferenceGenomeDifference.create(20, "G", "C")),
              Arrays.asList(),
              Arrays.asList(
                  TransformationUnit.create(
                      ReferenceGenomeDifference.create(17, "A", "G"),
                      Arrays.<Variant>asList(),
                      LARGEST_CHR1_HOMOLOG),
                  TransformationUnit.create(
                      ReferenceGenomeDifference.create(20, "G", "C"),
                      Arrays.<Variant>asList(),
                      LARGEST_CHR1_HOMOLOG))
            },
            // Both reference genome changes and variants.
            {
              Arrays.asList(ReferenceGenomeDifference.create(17, "A", "G")),
              Arrays.asList(SINGLE_NUCLEOTIDE_HET1),
              Arrays.asList(
                  TransformationUnit.create(
                      ReferenceGenomeDifference.create(17, "A", "G"),
                      Arrays.asList(SINGLE_NUCLEOTIDE_HET1),
                      LARGEST_CHR1_HOMOLOG))
            },
            {
              Arrays.asList(ReferenceGenomeDifference.create(20, "G", "C")),
              Arrays.asList(SINGLE_NUCLEOTIDE_HET1),
              Arrays.asList(
                  TransformationUnit.create(
                      ReferenceGenomeDifference.create(20, "G", "C"),
                      Arrays.<Variant>asList(),
                      LARGEST_CHR1_HOMOLOG),
                  TransformationUnit.create(
                      ReferenceGenomeDifference.NO_DIFFERENCE,
                      Arrays.asList(SINGLE_NUCLEOTIDE_HET1),
                      LARGEST_CHR1_HOMOLOG))
            }
          });
    }

    @Parameter(0)
    public List<ReferenceGenomeDifference> referenceGenomeDifferences;

    @Parameter(1)
    public List<Variant> variants;

    @Parameter(2)
    public List<TransformationUnit> units;

    @Test
    public void createTransformationUnits() {
      List<TransformationUnit> actual =
          TransformVariantsInRange.createTransformationUnits(
              referenceGenomeDifferences, variants, LARGEST_CHR1_HOMOLOG);
      assertTrue(actual.size() == units.size());
      List<TransformationUnit> mutableUnits = new ArrayList<>(units);
      for (TransformationUnit actualUnit : actual) {
        for (int i = 0; i < mutableUnits.size(); i++) {
          if (actualUnit
              .referenceGenomeDifference()
              .equals(mutableUnits.get(i).referenceGenomeDifference())) {
            TransformationUnit expectedUnit = mutableUnits.remove(i);
            List<Variant> actualVariants = actualUnit.variants();
            List<Variant> expectedVariants = expectedUnit.variants();
            assertTrue(actualVariants.containsAll(expectedVariants)
                && expectedVariants.containsAll(actualVariants));
          }
        }
      }
      assertTrue(mutableUnits.size() == 0);
    }
  }

  public static final class CreateUnsupportedTransformationUnitsTests {

    @Test
    public void createInvalidTransformationUnits() {
      List<ReferenceGenomeDifference> differences =
          Arrays.asList(
              ReferenceGenomeDifference.create(27, "C", "T"),
              ReferenceGenomeDifference.create(29, "T", "A"));
      List<Variant> overlappingVariants =
          Arrays.asList(
              makeVariant("chr1", 27, "CATG", "C", createSomeCalls(CALLSET_NAME, 1, 1)),
              makeVariant("chr1", 29, "T", "C", createSomeCalls(CALLSET_NAME, 0, 0)));

      List<TransformationUnit> actual =
          TransformVariantsInRange.createTransformationUnits(
              differences, overlappingVariants, LARGEST_CHR1_HOMOLOG);
      assertTrue(actual == null); // The sentinel value for INVALID_TRANSFORMATION_UNITS.
    }
  }

  @RunWith(Parameterized.class)
  public static final class TargetVariantsForAllTransformationUnitsTests {
    @Parameters(name = "{index}: targetVariantsForAllTransformationUnits({0})={1}")
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            {Arrays.asList(), Arrays.asList()},
            {
              Arrays.asList(
                  TransformationUnit.create(
                      ReferenceGenomeDifference.NO_DIFFERENCE,
                      Arrays.asList(SINGLE_NUCLEOTIDE_HET1),
                      LARGEST_CHR1_HOMOLOG)),
              Arrays.asList(makeVariant("chr1_same", 27, "A", "G",
                  createSomeCalls(CALLSET_NAME, 0, 1, 1, 1, 1, 0)))
            },
            // Checks the retention of homozygous reference target variants.
            {
              Arrays.asList(
                  TransformationUnit.create(
                      ReferenceGenomeDifference.create(3, "G", "A"),
                      Arrays.asList(makeVariant("chr1", 3, "G", "A",
                        createSomeCalls(CALLSET_NAME, 1, 1))),
                      LARGEST_CHR1_HOMOLOG)),
              Arrays.asList(makeVariant("chr1_same", 13, "A", "G",
                  createSomeCalls(CALLSET_NAME, 0, 0)))
            },
            {
              Arrays.asList(
                  TransformationUnit.create(
                      ReferenceGenomeDifference.create(3, "G", "A"),
                      Arrays.asList(makeVariant("chr1", 3, "G", "A",
                        createSomeCalls(CALLSET_NAME, 1, 1, 1, 0))),
                      LARGEST_CHR1_HOMOLOG),
                  TransformationUnit.create(
                      ReferenceGenomeDifference.NO_DIFFERENCE,
                      Arrays.asList(SINGLE_NUCLEOTIDE_HET1),
                      LARGEST_CHR1_HOMOLOG)),
              Arrays.asList(
                  makeVariant("chr1_same", 13, "A", "G", createSomeCalls(CALLSET_NAME, 0, 0, 0, 1)),
                  makeVariant("chr1_same", 27, "A", "G",
                    createSomeCalls(CALLSET_NAME, 0, 1, 1, 1, 1, 0)))
            }
          });
    }

    @Parameter(0)
    public List<TransformationUnit> units;

    @Parameter(1)
    public List<Variant> expectedVariants;

    @Test
    public void validTransformationUnitsTests() {
      TransformedVariants targetVariants =
          TransformVariantsInRange.getTargetVariantsForAllTransformationUnitsInRange(
              units, CALLSET_NAME);
      assertTrue(targetVariants.getStatus().equals(TransformationStatus.SUCCESS));
      List<Variant> actualVariants = targetVariants.getVariantsList();
      assertTrue(actualVariants.containsAll(expectedVariants)
          && expectedVariants.containsAll(actualVariants));
    }
  }

  @RunWith(Parameterized.class)
  public static final class UnsupportedTargetVariantsForAllTransformationUnitsTests {

    private static final TransformationUnit UNSUPPORTED_CONVERSION_UNIT =
        TransformationUnit.create(
            ReferenceGenomeDifference.NO_DIFFERENCE,
            Arrays.asList(INSERTION_HET1),
            LARGEST_CHR2_REVCOMP_HOMOLOG);

    @Parameters(name = "{index}: unsupportedGetTargetVariantsForAllTransformationUnits({0})")
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            {null}, // The sentinel object for invalid transformation units.
            {Arrays.asList(UNSUPPORTED_CONVERSION_UNIT)},
            {
              Arrays.asList(
                  TransformationUnit.create(
                      ReferenceGenomeDifference.NO_DIFFERENCE,
                      Arrays.asList(SINGLE_NUCLEOTIDE_HET1),
                      LARGEST_CHR1_HOMOLOG),
                  UNSUPPORTED_CONVERSION_UNIT)
            }
          });
    }

    @Parameter(0)
    public List<TransformationUnit> units;

    @Test
    public void validTransformationUnitsTests() {
      TransformedVariants targetVariants =
          TransformVariantsInRange.getTargetVariantsForAllTransformationUnitsInRange(
              units, CALLSET_NAME);
      assertTrue(targetVariants.getStatus().equals(TransformationStatus.FAILURE));
      assertTrue(targetVariants.getVariantsList().size() == 0);
    }
  }

  @RunWith(Parameterized.class)
  public static final class ReferenceChangesFromIndelTests {
    @Parameters(name = "{index}: posStrandIndelChange({0}, {1})={2}")
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            // Indels that are not replicating the exact sequence lost don't have differences.
            {
              LARGEST_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 1, "AC", "A", createSomeCalls(CALLSET_NAME, 0, 1, 1, 1)),
              ReferenceGenomeDifference.NO_DIFFERENCE
            },
            {
              LARGEST_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 1, "A", "ACT", createSomeCalls(CALLSET_NAME, 0, 1)),
              ReferenceGenomeDifference.NO_DIFFERENCE
            },
            {
              LARGEST_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 1, "ACTGC", "A", createSomeCalls(CALLSET_NAME, 0, 1)),
              ReferenceGenomeDifference.NO_DIFFERENCE
            },
            {
              LARGEST_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 1, "ACTGCT", "A", createSomeCalls(CALLSET_NAME, 0, 1)),
              ReferenceGenomeDifference.NO_DIFFERENCE
            },
            {
              LARGEST_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 1, "ACTGCTG", "A", createSomeCalls(CALLSET_NAME, 0, 1)),
              ReferenceGenomeDifference.NO_DIFFERENCE
            },
            // The exact indel detects a copy number change.
            {
              LARGEST_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 1, "ACTG", "A", createSomeCalls(CALLSET_NAME, 0, 1)),
              ReferenceGenomeDifference.create(1, "ACTG", "A")
            },
            {
              LARGEST_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 1, "A", "ACTG", createSomeCalls(CALLSET_NAME, 0, 1)),
              ReferenceGenomeDifference.create(1, "ACTG", "A")
            },
            // All the above tests should also work for a shorter homologous region in which we peek
            // out of it.
            {
              SMALL_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 1, "AC", "A", createSomeCalls(CALLSET_NAME, 0, 1)),
              ReferenceGenomeDifference.NO_DIFFERENCE
            },
            {
              SMALL_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 1, "A", "ACT", createSomeCalls(CALLSET_NAME, 0, 1)),
              ReferenceGenomeDifference.NO_DIFFERENCE
            },
            {
              SMALL_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 1, "ACTGC", "A", createSomeCalls(CALLSET_NAME, 0, 1)),
              ReferenceGenomeDifference.NO_DIFFERENCE
            },
            {
              SMALL_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 1, "ACTGCT", "A", createSomeCalls(CALLSET_NAME, 0, 1)),
              ReferenceGenomeDifference.NO_DIFFERENCE
            },
            {
              SMALL_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 1, "ACTGCTG", "A", createSomeCalls(CALLSET_NAME, 0, 1)),
              ReferenceGenomeDifference.NO_DIFFERENCE
            },
            {
              SMALL_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 1, "ACTG", "A", createSomeCalls(CALLSET_NAME, 0, 1)),
              ReferenceGenomeDifference.create(1, "ACTG", "A")
            },
            {
              SMALL_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 1, "A", "ACTG", createSomeCalls(CALLSET_NAME, 0, 1)),
              ReferenceGenomeDifference.create(1, "ACTG", "A")
            },
            // Multiallelic variants containing the indel should also detect a change.
            {
              SMALL_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              Variant.newBuilder()
                  .setReferenceName("chr2")
                  .setStart(1)
                  .setEnd(2)
                  .setReferenceBases("A")
                  .addAlternateBases("AC")
                  .addAlternateBases("ACTG")
                  .addAllCalls(createSomeCalls(CALLSET_NAME, 0, 1))
                  .build(),
              ReferenceGenomeDifference.create(1, "ACTG", "A")
            },
            {
              LARGEST_CHR2_CTG_INSERTION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 1, "A", "ACTG", createSomeCalls(CALLSET_NAME, 0, 1)),
              ReferenceGenomeDifference.create(1, "A", "ACTG")
            },
            {
              LARGEST_CHR2_CTG_INSERTION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 21, "T", "TGCT", createSomeCalls(CALLSET_NAME, 0, 1)),
              ReferenceGenomeDifference.create(21, "T", "TGCT")
            },
            {
              LARGEST_CHR2_CTG_INSERTION_IDENTICAL_HOMOLOG,
              makeVariant("chr2", 2, "TGC", "GGGG", createSomeCalls(CALLSET_NAME, 0, 1)),
              ReferenceGenomeDifference.NO_DIFFERENCE
            }
          });
    }

    @Parameter(0)
    public HomologousRange range;

    @Parameter(1)
    public Variant indel;

    @Parameter(2)
    public ReferenceGenomeDifference expected;

    @Test
    public void posStrandIndelChange() throws IOException {
      Fasta queryFasta = new Fasta(TransformVariantsInRangeTest.class.getClassLoader()
          .getResource(QUERY_REFERENCE_PATH).getFile());
      Fasta targetFasta = new Fasta(TransformVariantsInRangeTest.class.getClassLoader()
          .getResource(TARGET_REFERENCE_PATH).getFile());
      ReferenceGenomeDifference actual =
          TransformVariantsInRange.getReferenceChangeFromIndelInPositiveStrandIdenticalRegion(
              range, indel, queryFasta, targetFasta);
      assertTrue(actual.equals(expected));
    }
  }
}
