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
import com.verily.genomewarp.TransformationUnit;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.RegionType;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.TargetStrand;
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
public final class TransformationUnitTest {

  private static final List<String> CALLSET_NAME = Arrays.asList("MyCallsetName1",
      "MyCallsetName2", "MyCallsetName3");

  private static final HomologousRange POSITIVE_STRAND_IDENTICAL_RANGE =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr1").setStart(12).setEnd(34))
          .setTargetRange(Range.newBuilder().setReferenceName("chr10").setStart(42).setEnd(64))
          .setTargetStrand(TargetStrand.POSITIVE_STRAND)
          .setRegionType(RegionType.IDENTICAL)
          .build();

  private static final HomologousRange NEGATIVE_STRAND_IDENTICAL_RANGE =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr1").setStart(12).setEnd(34))
          .setTargetRange(Range.newBuilder().setReferenceName("chr10").setStart(42).setEnd(64))
          .setTargetStrand(TargetStrand.NEGATIVE_STRAND)
          .setRegionType(RegionType.IDENTICAL)
          .build();

  private static final HomologousRange POSITIVE_STRAND_MISMATCH_RANGE =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr1").setStart(12).setEnd(34))
          .setTargetRange(Range.newBuilder().setReferenceName("chr10").setStart(42).setEnd(64))
          .setTargetStrand(TargetStrand.POSITIVE_STRAND)
          .setRegionType(RegionType.MISMATCHED_BASES)
          .build();

  private static final HomologousRange NEGATIVE_STRAND_MISMATCH_RANGE =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr1").setStart(12).setEnd(34))
          .setTargetRange(Range.newBuilder().setReferenceName("chr10").setStart(42).setEnd(64))
          .setTargetStrand(TargetStrand.NEGATIVE_STRAND)
          .setRegionType(RegionType.MISMATCHED_BASES)
          .build();

  private static final HomologousRange NO_STRAND_RANGE =
      HomologousRange.newBuilder()
          .setQueryRange(Range.newBuilder().setReferenceName("chr1").setStart(12).setEnd(34))
          .setTargetRange(Range.newBuilder().setReferenceName("chr10").setStart(42).setEnd(64))
          .build();

  private static final Variant SINGLE_NUCLEOTIDE_HOMREF =
      makeVariant("chr1", 17, "A", "G", createSomeCalls(CALLSET_NAME, 0, 0, 1, 1));
  private static final Variant SINGLE_NUCLEOTIDE_HET1 =
      makeVariant("chr1", 17, "A", "G", createSomeCalls(CALLSET_NAME, 0, 1, 1, 0, 1, 1));
  private static final Variant SINGLE_NUCLEOTIDE_HET2 =
      makeVariant("chr1", 17, "A", "G", createSomeCalls(CALLSET_NAME, 1, 0));
  private static final Variant SINGLE_NUCLEOTIDE_HOMALT =
      makeVariant("chr1", 17, "A", "G", createSomeCalls(CALLSET_NAME, 1, 1));

  private static final Variant INSERTION_HOMREF =
      makeVariant("chr1", 20, "A", "AT", createSomeCalls(CALLSET_NAME, 0, 0));
  private static final Variant INSERTION_HET1 =
      makeVariant("chr1", 20, "A", "AT", createSomeCalls(CALLSET_NAME, 0, 1));
  private static final Variant INSERTION_HET2 =
      makeVariant("chr1", 20, "A", "AT", createSomeCalls(CALLSET_NAME, 1, 0));
  private static final Variant INSERTION_HOMALT =
      makeVariant("chr1", 20, "A", "AT", createSomeCalls(CALLSET_NAME, 1, 1));

  private static final Variant DELETION_HOMREF =
      makeVariant("chr1", 23, "CC", "C", createSomeCalls(CALLSET_NAME, 0, 0, 1, 0));
  private static final Variant DELETION_HET1 =
      makeVariant("chr1", 23, "CC", "C", createSomeCalls(CALLSET_NAME, 0, 1));
  private static final Variant DELETION_HET2 =
      makeVariant("chr1", 23, "CC", "C", createSomeCalls(CALLSET_NAME, 1, 0));
  private static final Variant DELETION_HOMALT =
      makeVariant("chr1", 23, "CC", "C", createSomeCalls(CALLSET_NAME, 1, 1));

  @RunWith(Parameterized.class)
  public static final class ValidCreateTests {
    @Parameters(name = "{index}: create({0})={1}")
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            {
              ReferenceGenomeDifference.NO_DIFFERENCE,
              Arrays.asList(),
              POSITIVE_STRAND_IDENTICAL_RANGE
            },
            {
              ReferenceGenomeDifference.NO_DIFFERENCE,
              Arrays.asList(SINGLE_NUCLEOTIDE_HOMREF),
              POSITIVE_STRAND_IDENTICAL_RANGE
            },
            {
              ReferenceGenomeDifference.NO_DIFFERENCE,
              Arrays.asList(INSERTION_HET1),
              NEGATIVE_STRAND_IDENTICAL_RANGE
            },
            {
              ReferenceGenomeDifference.NO_DIFFERENCE,
              Arrays.asList(DELETION_HET2),
              NEGATIVE_STRAND_IDENTICAL_RANGE
            },
            {
              ReferenceGenomeDifference.NO_DIFFERENCE,
              Arrays.asList(DELETION_HOMALT, DELETION_HOMREF),
              POSITIVE_STRAND_IDENTICAL_RANGE
            },
            {
              ReferenceGenomeDifference.create(15, "A", "G"),
              Arrays.asList(),
              POSITIVE_STRAND_MISMATCH_RANGE
            },
            {
              ReferenceGenomeDifference.create(17, "A", "G"),
              Arrays.asList(SINGLE_NUCLEOTIDE_HOMALT),
              NEGATIVE_STRAND_MISMATCH_RANGE
            },
          });
    }

    @Parameter(0)
    public ReferenceGenomeDifference expectedReferenceGenomeDifference;

    @Parameter(1)
    public List<Variant> expectedVariants;

    @Parameter(2)
    public HomologousRange expectedRange;

    @Test
    public void createTransformationUnit() {
      TransformationUnit actualUnit =
          TransformationUnit.create(
              expectedReferenceGenomeDifference, expectedVariants, expectedRange);
      assertTrue(actualUnit.referenceGenomeDifference().equals(expectedReferenceGenomeDifference));
      assertTrue(actualUnit.variants().equals(expectedVariants));
      assertTrue(actualUnit.range().equals(expectedRange));
    }
  }

  @RunWith(Parameterized.class)
  public static final class SimpleCoordinateTransformTests {

    @Parameters(name = "{index}: simpleCoordinateTransform({0}, {1})={2}")
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            // SNVs on the positive strand.
            {
              POSITIVE_STRAND_IDENTICAL_RANGE,
              SINGLE_NUCLEOTIDE_HOMREF,
              makeVariant("chr10", 47, "A", "G", createSomeCalls(CALLSET_NAME, 0, 0, 1, 1))
            },
            {
              POSITIVE_STRAND_IDENTICAL_RANGE,
              SINGLE_NUCLEOTIDE_HET1,
              makeVariant("chr10", 47, "A", "G", createSomeCalls(CALLSET_NAME, 0, 1, 1, 0, 1, 1))
            },
            {
              POSITIVE_STRAND_IDENTICAL_RANGE,
              SINGLE_NUCLEOTIDE_HET2,
              makeVariant("chr10", 47, "A", "G", createSomeCalls(CALLSET_NAME, 1, 0))
            },
            {
              POSITIVE_STRAND_IDENTICAL_RANGE,
              SINGLE_NUCLEOTIDE_HOMALT,
              makeVariant("chr10", 47, "A", "G", createSomeCalls(CALLSET_NAME, 1, 1))
            },
            // Insertions on the positive strand.
            {
              POSITIVE_STRAND_IDENTICAL_RANGE,
              INSERTION_HOMREF,
              makeVariant("chr10", 50, "A", "AT", createSomeCalls(CALLSET_NAME, 0, 0))
            },
            {
              POSITIVE_STRAND_IDENTICAL_RANGE,
              INSERTION_HET1,
              makeVariant("chr10", 50, "A", "AT", createSomeCalls(CALLSET_NAME, 0, 1))
            },
            {
              POSITIVE_STRAND_IDENTICAL_RANGE,
              INSERTION_HET2,
              makeVariant("chr10", 50, "A", "AT", createSomeCalls(CALLSET_NAME, 1, 0))
            },
            {
              POSITIVE_STRAND_IDENTICAL_RANGE,
              INSERTION_HOMALT,
              makeVariant("chr10", 50, "A", "AT", createSomeCalls(CALLSET_NAME, 1, 1))
            },
            // Deletions on the positive strand.
            {
              POSITIVE_STRAND_IDENTICAL_RANGE,
              DELETION_HOMREF,
              makeVariant("chr10", 53, "CC", "C", createSomeCalls(CALLSET_NAME, 0, 0, 1, 0))
            },
            {
              POSITIVE_STRAND_IDENTICAL_RANGE,
              DELETION_HET1,
              makeVariant("chr10", 53, "CC", "C", createSomeCalls(CALLSET_NAME, 0, 1))
            },
            {
              POSITIVE_STRAND_IDENTICAL_RANGE,
              DELETION_HET2,
              makeVariant("chr10", 53, "CC", "C", createSomeCalls(CALLSET_NAME, 1, 0))
            },
            {
              POSITIVE_STRAND_IDENTICAL_RANGE,
              DELETION_HOMALT,
              makeVariant("chr10", 53, "CC", "C", createSomeCalls(CALLSET_NAME, 1, 1))
            },
            // SNVs on the negative strand.
            {
              NEGATIVE_STRAND_IDENTICAL_RANGE,
              SINGLE_NUCLEOTIDE_HOMREF,
              makeVariant("chr10", 58, "T", "C", createSomeCalls(CALLSET_NAME, 0, 0, 1, 1))
            },
            {
              NEGATIVE_STRAND_IDENTICAL_RANGE,
              SINGLE_NUCLEOTIDE_HET1,
              makeVariant("chr10", 58, "T", "C", createSomeCalls(CALLSET_NAME, 0, 1, 1, 0, 1, 1))
            },
            {
              NEGATIVE_STRAND_IDENTICAL_RANGE,
              SINGLE_NUCLEOTIDE_HET2,
              makeVariant("chr10", 58, "T", "C", createSomeCalls(CALLSET_NAME, 1, 0))
            },
            {
              NEGATIVE_STRAND_IDENTICAL_RANGE,
              SINGLE_NUCLEOTIDE_HOMALT,
              makeVariant("chr10", 58, "T", "C", createSomeCalls(CALLSET_NAME, 1, 1))
            },
          });
    }

    @Parameter(0)
    public HomologousRange range;

    @Parameter(1)
    public Variant queryVariant;

    @Parameter(2)
    public Variant expectedTargetVariant;

    @Test
    public void performSimpleCoordinateTransform() {
      Variant actualTargetVariant =
          TransformationUnit.simpleCoordinateTransform(
              range,
              queryVariant,
              queryVariant.getStart(),
              queryVariant.getReferenceBases(),
              queryVariant.getAlternateBasesList());
      assertTrue(actualTargetVariant.equals(expectedTargetVariant));
    }
  }

  public static final class InvalidSimpleCoordinateTransformTests {

    @Rule public final ExpectedException thrown = ExpectedException.none();

    @Test
    public void invalidSimpleCoordinateTransform() {
      thrown.expect(IllegalArgumentException.class);
      TransformationUnit.simpleCoordinateTransform(
          NO_STRAND_RANGE,
          makeVariant("chr10", 58, "T", "C", createSomeCalls(CALLSET_NAME, 1, 1)),
          12,
          "T",
          Arrays.asList("C"));
    }
  }

  @RunWith(Parameterized.class)
  public static final class GetTargetAssemblyVariantsTests {

    @Parameters(name = "{index}: {0}.getTargetAssemblyVariants()={1}")
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            // Variants in regions with no reference genome differences.
            {
              TransformationUnit.create(
                  ReferenceGenomeDifference.NO_DIFFERENCE,
                  Arrays.<Variant>asList(),
                  POSITIVE_STRAND_IDENTICAL_RANGE),
              Arrays.asList()
            },
            {
              TransformationUnit.create(
                  ReferenceGenomeDifference.NO_DIFFERENCE,
                  Arrays.asList(SINGLE_NUCLEOTIDE_HOMREF),
                  POSITIVE_STRAND_IDENTICAL_RANGE),
              Arrays.asList(makeVariant("chr10", 47, "A", "G",
                  createSomeCalls(CALLSET_NAME, 0, 0, 1, 1)))
            },
            {
              TransformationUnit.create(
                  ReferenceGenomeDifference.NO_DIFFERENCE,
                  Arrays.asList(SINGLE_NUCLEOTIDE_HOMREF, INSERTION_HET1),
                  POSITIVE_STRAND_IDENTICAL_RANGE),
              Arrays.asList(
                  makeVariant("chr10", 47, "A", "G", createSomeCalls(CALLSET_NAME, 0, 0, 1, 1)),
                  makeVariant("chr10", 50, "A", "AT", createSomeCalls(CALLSET_NAME, 0, 1)))
            },
            {
              TransformationUnit.create(
                  ReferenceGenomeDifference.NO_DIFFERENCE,
                  Arrays.asList(SINGLE_NUCLEOTIDE_HOMREF, SINGLE_NUCLEOTIDE_HET2),
                  NEGATIVE_STRAND_IDENTICAL_RANGE),
              Arrays.asList(
                  makeVariant("chr10", 58, "T", "C", createSomeCalls(CALLSET_NAME, 0, 0, 1, 1)),
                  makeVariant("chr10", 58, "T", "C", createSomeCalls(CALLSET_NAME, 1, 0)))
            },
            // Reference-genome-only differences
            {
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(20, "A", "C"),
                  Arrays.<Variant>asList(),
                  POSITIVE_STRAND_MISMATCH_RANGE),
              Arrays.asList(makeVariant("chr10", 50, "C", "A", createSomeCalls(CALLSET_NAME, 1, 1)))
            },
            {
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(21, "A", "C"),
                  Arrays.<Variant>asList(),
                  NEGATIVE_STRAND_MISMATCH_RANGE),
              Arrays.asList(makeVariant("chr10", 54, "G", "T", createSomeCalls(CALLSET_NAME, 1, 1)))
            },
            {
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(20, "A", "AC"),
                  Arrays.<Variant>asList(),
                  POSITIVE_STRAND_MISMATCH_RANGE),
              Arrays.asList(makeVariant("chr10", 50, "AC", "A",
                  createSomeCalls(CALLSET_NAME, 1, 1)))
            },
            // Single nucleotide changes in both genome reference and the variant
            {
              // Reference changes from A --> G, on the A reference we are homozygous GG
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(17, "A", "G"),
                  Arrays.<Variant>asList(SINGLE_NUCLEOTIDE_HOMALT),
                  POSITIVE_STRAND_MISMATCH_RANGE),
              Arrays.asList(makeVariant("chr10", 47, "G", "A", createSomeCalls(CALLSET_NAME, 0, 0)))
            },
            {
              // Reference changes from A --> G, on the A reference we are heterozygous G/A
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(17, "A", "G"),
                  Arrays.<Variant>asList(SINGLE_NUCLEOTIDE_HET2),
                  POSITIVE_STRAND_MISMATCH_RANGE),
              Arrays.asList(
                  makeVariant(
                      "chr10",
                      47,
                      "G",
                      "A",
                      // Note that the ordering preserves "G" in the first position
                      createSomeCalls(CALLSET_NAME, 0, 1)))
            },
            {
              // Reference changes from A --> G, on the A reference is an A/G variant but it is a
              // no-call
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(17, "A", "G"),
                  Arrays.asList(makeVariant("chr1", 17, "A", "G",
                    createSomeCalls(CALLSET_NAME, -1, -1))),
                  POSITIVE_STRAND_MISMATCH_RANGE),
              Arrays.asList(makeVariant("chr10", 47, "G", "A",
                  createSomeCalls(CALLSET_NAME, -1, -1)))
            },
            {
              // Reference changes from A --> T, on the A reference we are heterozygous A/G
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(17, "A", "T"),
                  Arrays.<Variant>asList(SINGLE_NUCLEOTIDE_HET1),
                  POSITIVE_STRAND_MISMATCH_RANGE),
              Arrays.asList(
                  Variant.newBuilder()
                      .setReferenceName("chr10")
                      .setStart(47)
                      .setEnd(48)
                      .setReferenceBases("T")
                      .addAlternateBases("A")
                      .addAlternateBases("G")
                      .addFilter("PASS")
                      .addAllCalls(createSomeCalls(CALLSET_NAME, 1, 2, 2, 1, 2, 2))
                      .build())
            },
            {
              // Reference changes from A --> G, on the A reference we are homozygous GG
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(17, "A", "G"),
                  Arrays.<Variant>asList(SINGLE_NUCLEOTIDE_HOMALT),
                  NEGATIVE_STRAND_MISMATCH_RANGE),
              Arrays.asList(makeVariant("chr10", 58, "C", "T", createSomeCalls(CALLSET_NAME, 0, 0)))
            },
            {
              // Reference changes from A --> G, on the A reference we are heterozygous A/G
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(17, "A", "G"),
                  Arrays.<Variant>asList(SINGLE_NUCLEOTIDE_HET1),
                  NEGATIVE_STRAND_MISMATCH_RANGE),
              Arrays.asList(makeVariant("chr10", 58, "C", "T",
                  createSomeCalls(CALLSET_NAME, 1, 0, 0, 1, 0, 0)))
            },
            {
              // Reference changes from A --> C, on the A reference we are heterozygous G/A
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(17, "A", "C"),
                  Arrays.<Variant>asList(SINGLE_NUCLEOTIDE_HET2),
                  NEGATIVE_STRAND_MISMATCH_RANGE),
              Arrays.asList(
                  Variant.newBuilder()
                      .setReferenceName("chr10")
                      .setStart(58)
                      .setEnd(59)
                      .setReferenceBases("G")
                      .addAlternateBases("T")
                      .addAlternateBases("C")
                      .addFilter("PASS")
                      .addAllCalls(createSomeCalls(CALLSET_NAME, 2, 1))
                      .build())
            },
            // Insertion/deletion variant tests
            {
              // Reference changes from A --> AT, on the A reference we are het AT/A
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(20, "A", "AT"),
                  Arrays.<Variant>asList(INSERTION_HET2),
                  POSITIVE_STRAND_IDENTICAL_RANGE),
              Arrays.asList(makeVariant("chr10", 50, "AT", "A",
                  createSomeCalls(CALLSET_NAME, 0, 1)))
            },
            {
              // Reference changes from A --> AT, on the A reference is an A/AT variant but this is
              // a no-call
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(20, "A", "AT"),
                  Arrays.asList(makeVariant("chr1", 20, "A", "AT",
                    createSomeCalls(CALLSET_NAME, -1, -1))),
                  POSITIVE_STRAND_IDENTICAL_RANGE),
              Arrays.asList(makeVariant("chr10", 50, "AT", "A",
                  createSomeCalls(CALLSET_NAME, -1, -1)))
            },
          });
    }

    @Parameter(0)
    public TransformationUnit unit;

    @Parameter(1)
    public List<Variant> expectedTargetVariants;

    @Test
    public void getTargetAssemblyVariants() {
      List<Variant> actualTargetVariants = unit.getTargetAssemblyVariants(CALLSET_NAME);
      assertEquals(actualTargetVariants, expectedTargetVariants);
    }
  }

  @RunWith(Parameterized.class)
  public static final class UnsupportedGetTargetAssemblyVariantsTests {

    @Parameters(name = "{index}: unsupportedGetTargetAssemblyVariants({0})=INVALID_TARGET_VARIANTS")
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            // Alignment-required region
            {
              TransformationUnit.create(
                  ReferenceGenomeDifference.NO_DIFFERENCE,
                  Arrays.asList(INSERTION_HET1),
                  HomologousRange.newBuilder()
                      .setQueryRange(
                          Range.newBuilder().setReferenceName("chr1").setStart(12).setEnd(34))
                      .setTargetRange(
                          Range.newBuilder().setReferenceName("chr10").setStart(42).setEnd(63))
                      .setTargetStrand(TargetStrand.POSITIVE_STRAND)
                      .setRegionType(RegionType.ALIGNMENT_REQUIRED)
                      .build())
            },
            // Insertion on the negative strand
            {
              TransformationUnit.create(
                  ReferenceGenomeDifference.NO_DIFFERENCE,
                  Arrays.asList(INSERTION_HET1),
                  NEGATIVE_STRAND_IDENTICAL_RANGE)
            },
            // Multiple variants for a single ReferenceGenomeDifference
            {
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(17, "A", "C"),
                  Arrays.<Variant>asList(SINGLE_NUCLEOTIDE_HOMREF, SINGLE_NUCLEOTIDE_HET2),
                  POSITIVE_STRAND_MISMATCH_RANGE)
            },
            {
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(17, "A", "C"),
                  Arrays.<Variant>asList(SINGLE_NUCLEOTIDE_HOMREF, SINGLE_NUCLEOTIDE_HET1),
                  NEGATIVE_STRAND_MISMATCH_RANGE)
            },
            // Deletion on the negative strand
            {
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(23, "A", "AT"),
                  Arrays.<Variant>asList(DELETION_HET2),
                  NEGATIVE_STRAND_IDENTICAL_RANGE)
            },
            // Mismatching indel between reference and variant
            {
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(20, "A", "AC"),
                  Arrays.<Variant>asList(INSERTION_HET2),
                  POSITIVE_STRAND_IDENTICAL_RANGE)
            },
            // Multiple alleles for variant
            {
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(20, "A", "AT"),
                  Arrays.<Variant>asList(
                      INSERTION_HET2.toBuilder().addAlternateBases("ATT").build()),
                  POSITIVE_STRAND_IDENTICAL_RANGE)
            },
            // Reference indel on the negative strand
            {
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(20, "AC", "A"),
                  Arrays.<Variant>asList(),
                  NEGATIVE_STRAND_IDENTICAL_RANGE)
            },
          });
    }

    @Parameter(0)
    public TransformationUnit unit;

    @Test
    public void unsupportedGetTargetAssemblyVariants() {
      List<Variant> actualTargetVariants = unit.getTargetAssemblyVariants(CALLSET_NAME);

      // Essentially ensures this is null
      assertTrue(actualTargetVariants == TransformationUnit.INVALID_TARGET_VARIANTS);
    }
  }

  @RunWith(Parameterized.class)
  public static final class InvalidGetTargetAssemblyVariantsTests {

    @Parameters(name = "{index}: invalidGetTargetAssemblyVariants({0})")
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            {
              // Reference changes from C --> G, but variant indicates A as reference
              TransformationUnit.create(
                  ReferenceGenomeDifference.create(17, "C", "G"),
                  Arrays.<Variant>asList(SINGLE_NUCLEOTIDE_HET2),
                  POSITIVE_STRAND_MISMATCH_RANGE)
            },
          });
    }

    @Parameter(0)
    public TransformationUnit unit;

    @Rule public final ExpectedException thrown = ExpectedException.none();

    @Test
    public void invalidTransformCluster() {
      thrown.expect(IllegalArgumentException.class);
      unit.getTargetAssemblyVariants(CALLSET_NAME);
    }
  }
}
