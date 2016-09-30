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
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.genomics.v1.Range;
import com.google.genomics.v1.Variant;
import com.google.genomics.v1.VariantCall;
import com.verily.genomewarp.GenomeWarpSerial;
import com.verily.genomewarp.TransformVariantsInRange;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.RegionType;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.TargetStrand;
import com.verily.genomewarp.utils.Fasta;
import com.verily.genomewarp.utils.GenomeRange;
import com.verily.genomewarp.utils.GenomeWarpTestUtils;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.junit.Rule;
import org.junit.Test;
import org.junit.experimental.runners.Enclosed;
import org.junit.rules.ExpectedException;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

@RunWith(Enclosed.class)
public final class GenomeWarpSerialTest {

  private static HomologousRange newHRange(Range a, Range b, boolean s, RegionType r) {
    TargetStrand t = s ? TargetStrand.POSITIVE_STRAND : TargetStrand.NEGATIVE_STRAND;
    HomologousRange.Builder builder = HomologousRange.newBuilder().setQueryRange(a)
        .setTargetRange(b).setTargetStrand(t);
    if (r != null) {
      builder.setRegionType(r);
    }

    return builder.build();
  }

  private static HomologousRange newHRange(Range a, Range b, boolean s) {
    return newHRange(a, b, s, null);
  }

  private static Range newRange(String ref, long start, long end) {
    return Range.newBuilder().setReferenceName(ref).setStart(start).setEnd(end).build();
  }

  @RunWith(JUnit4.class)
  public static final class UtilityTests {

    private static final String FASTA_PATH = "test.fa";


    private String fakeBED(String[][] in) {
      StringBuilder toReturn = new StringBuilder();
      for (int i = 0; i < in.length; i++) {
        for (int j = 0; j < in[i].length; j++) {
          toReturn.append(in[i][j] + "\t");
        }
        toReturn.append("\n");
      }
      return toReturn.toString();
    }

    @Test
    public void testSplitAtNonDNA() throws Exception {

      Fasta testFasta =
          new Fasta(GenomeWarpSerialTest.class.getClassLoader().getResource(FASTA_PATH).getFile());

      // Define mock BufferedReader
      String[][] inArray = {
        {"chr1", "10", "24"},
        {"chr2", "5", "6"},
        {"chr3", "9", "10"},
        {"chr4", "9", "19"},
        {"chr5", "10", "20"},
        {"chr6", "30", "40"},
        {"chr7", "40", "50"},
        {"chr8", "55", "72"},
        {"chr9", "10", "10"}};
      String testStringBED = fakeBED(inArray);
      InputStream stream = new ByteArrayInputStream(testStringBED.getBytes(StandardCharsets.UTF_8));
      BufferedReader testBED = new BufferedReader(new InputStreamReader(stream, UTF_8));

      // Define desired result
      GenomeRange[] wantArray = {
          new GenomeRange("chr1", 10, 24),
          new GenomeRange("chr2", 5, 6),
          new GenomeRange("chr5", 15, 20),
          new GenomeRange("chr6", 30, 35),
          new GenomeRange("chr7", 40, 44),
          new GenomeRange("chr7", 45, 50),
          new GenomeRange("chr8", 56, 60),
          new GenomeRange("chr8", 61, 66),
          new GenomeRange("chr8", 67, 72)};
      List<GenomeRange> want = Arrays.asList(wantArray);

      // Get test result
      List<GenomeRange> got =
          GenomeWarpSerial.splitAtNonDNA(testFasta, testBED);

      assertTrue(GenomeWarpTestUtils.equivalentRanges(got, want));
    }

    private void checkOmitOverlaps(GenomeRange[] in,
        GenomeRange[] want) {
      List<GenomeRange> inList, wantList;
      inList = Arrays.asList(in);
      wantList = Arrays.asList(want);

      List<GenomeRange> gotList =
          GenomeWarpSerial.omitOverlap(inList);

      assertTrue(GenomeWarpTestUtils.equivalentRanges(gotList, wantList));
    }

    @Test
    public void testOmitOverlapSorted() {
      GenomeRange[] in1 = {
          new GenomeRange("chr1", 10, 16),
          new GenomeRange("chr1", 20, 26),
          new GenomeRange("chr1", 30, 36),
          new GenomeRange("chr1", 40, 46),
          new GenomeRange("chr1", 50, 56),
          new GenomeRange("chr1", 60, 66),
          new GenomeRange("chr1", 70, 76),
          new GenomeRange("chr1", 80, 86),
          new GenomeRange("chr1", 90, 96)
      };
      GenomeRange[] want1 = {
          new GenomeRange("chr1", 10, 16),
          new GenomeRange("chr1", 20, 26),
          new GenomeRange("chr1", 30, 36),
          new GenomeRange("chr1", 40, 46),
          new GenomeRange("chr1", 50, 56),
          new GenomeRange("chr1", 60, 66),
          new GenomeRange("chr1", 70, 76),
          new GenomeRange("chr1", 80, 86),
          new GenomeRange("chr1", 90, 96)
      };
      checkOmitOverlaps(in1, want1);

      GenomeRange[] in2 = {
          new GenomeRange("chr1", 10, 16),
          new GenomeRange("chr1", 20, 26),
          new GenomeRange("chr1", 30, 79),
          new GenomeRange("chr1", 40, 46),
          new GenomeRange("chr1", 50, 56),
          new GenomeRange("chr1", 60, 66),
          new GenomeRange("chr1", 70, 76),
          new GenomeRange("chr1", 80, 86),
          new GenomeRange("chr1", 90, 96)
      };
      GenomeRange[] want2 = {
          new GenomeRange("chr1", 10, 16),
          new GenomeRange("chr1", 20, 26),
          new GenomeRange("chr1", 80, 86),
          new GenomeRange("chr1", 90, 96)
      };
      checkOmitOverlaps(in2, want2);

      GenomeRange[] in3 = {
          new GenomeRange("chr1", 10, 21),
          new GenomeRange("chr1", 20, 26),
          new GenomeRange("chr1", 30, 50),
          new GenomeRange("chr1", 40, 46),
          new GenomeRange("chr1", 50, 56),
          new GenomeRange("chr1", 60, 66),
          new GenomeRange("chr1", 70, 82),
          new GenomeRange("chr1", 80, 86),
          new GenomeRange("chr2", 10, 21),
          new GenomeRange("chr3", 10, 21),
          new GenomeRange("chr3", 100, 210),
          new GenomeRange("chr3", 220, 230),
          new GenomeRange("chr3", 225, 230)
      };
      GenomeRange[] want3 = {
          new GenomeRange("chr1", 50, 56),
          new GenomeRange("chr1", 60, 66),
          new GenomeRange("chr2", 10, 21),
          new GenomeRange("chr3", 10, 21),
          new GenomeRange("chr3", 100, 210)
      };
      checkOmitOverlaps(in3, want3);

      GenomeRange[] in4 = {
          new GenomeRange("chr1", 10, 21),
          new GenomeRange("chr1", 20, 26),
          new GenomeRange("chr1", 30, 51),
          new GenomeRange("chr1", 40, 46),
          new GenomeRange("chr1", 50, 56),
          new GenomeRange("chr1", 60, 66),
          new GenomeRange("chr1", 70, 82),
          new GenomeRange("chr1", 80, 86),
          new GenomeRange("chr2", 10, 21),
          new GenomeRange("chr2", 15, 17),
          new GenomeRange("chr3", 10, 21),
          new GenomeRange("chr3", 100, 210),
          new GenomeRange("chr3", 220, 230),
          new GenomeRange("chr3", 225, 230),
          new GenomeRange("chrX", 225, 230)
      };
      GenomeRange[] want4 = {
          new GenomeRange("chr1", 60, 66),
          new GenomeRange("chr3", 10, 21),
          new GenomeRange("chr3", 100, 210),
          new GenomeRange("chrX", 225, 230)
      };
      checkOmitOverlaps(in4, want4);
    }

    @Rule
    public ExpectedException thrown = ExpectedException.none();

    private void checkOmitOverlapsBad(GenomeRange[] in) {
      List<GenomeRange> inList = Arrays.asList(in);

      GenomeWarpSerial.omitOverlap(inList);
    }

    @Test
    public void testOmitOverlapUnsortedPos() {
      thrown.expect(RuntimeException.class);
      thrown.expectMessage("Output bed of liftover is not sorted by pos");
      GenomeRange[] in = {
          new GenomeRange("chr1", 10, 16),
          new GenomeRange("chr1", 20, 26),
          new GenomeRange("chr1", 40, 46),
          new GenomeRange("chr1", 30, 36),
          new GenomeRange("chr1", 50, 56),
          new GenomeRange("chr1", 60, 66),
          new GenomeRange("chr1", 70, 76),
          new GenomeRange("chr1", 80, 86),
          new GenomeRange("chr1", 90, 96)
      };
      checkOmitOverlapsBad(in);
    }

    @Test
    public void testOmitOverlapUnsortedChr() {
      thrown.expect(RuntimeException.class);
      thrown.expectMessage("Output bed of liftover is not sorted by chr");
      GenomeRange[] in = {
          new GenomeRange("chr1", 10, 16),
          new GenomeRange("chr1", 20, 26),
          new GenomeRange("chr2", 30, 36),
          new GenomeRange("chr2", 40, 46),
          new GenomeRange("chr2", 50, 56),
          new GenomeRange("chr2", 60, 66),
          new GenomeRange("chr1", 70, 76),
          new GenomeRange("chr1", 80, 86),
          new GenomeRange("chr1", 90, 96)
      };
      checkOmitOverlapsBad(in);
    }

    private void checkRegion(AnalyzeRegionsTest t) {
      HomologousRange range = HomologousRange.newBuilder()
          .setQueryRange(
              Range.newBuilder().setReferenceName(t.qChr).setStart(t.qStart).setEnd(t.qEnd).build())
          .setTargetRange(
              Range.newBuilder().setReferenceName(t.tChr).setStart(t.tStart).setEnd(t.tEnd).build())
          .setTargetStrand(t.strand ? TargetStrand.POSITIVE_STRAND : TargetStrand.NEGATIVE_STRAND)
          .build();

      RegionType out = GenomeWarpSerial.getRegionType(range, t.qFasta, t.tFasta);

      assertTrue(out == t.region);
    }

    private static class AnalyzeRegionsTest {
      public String name;
      public String qChr;
      public long qStart;
      public long qEnd;
      public String tChr;
      public long tStart;
      public long tEnd;
      public boolean strand;
      public Fasta qFasta;
      public Fasta tFasta;
      public RegionType region;

      public AnalyzeRegionsTest(String name, String qChr, long qStart, long qEnd, String tChr,
          long tStart, long tEnd, boolean strand, Fasta qFasta,
          Fasta tFasta, RegionType region) {
        this.name = name;
        this.qChr = qChr;
        this.qStart = qStart;
        this.qEnd = qEnd;
        this.tChr = tChr;
        this.tStart = tStart;
        this.tEnd = tEnd;
        this.strand = strand;
        this.qFasta = qFasta;
        this.tFasta = tFasta;
        this.region = region;
      }
    }

    @Test
    public void testAnalyzeRegions() throws IOException{
      String fastaPath =
          GenomeWarpSerialTest.class.getClassLoader().getResource(FASTA_PATH).getFile();
      AnalyzeRegionsTest[] in = {
        new AnalyzeRegionsTest("elt.1", "chr10a", 2, 10, "chr20a", 12, 20, true,
            new Fasta(fastaPath), new Fasta(fastaPath),
            RegionType.IDENTICAL),
        new AnalyzeRegionsTest("elt.2", "chr10b", 1, 9, "chr20b", 10, 19, true,
            new Fasta(fastaPath), new Fasta(fastaPath),
            RegionType.ALIGNMENT_REQUIRED),
        new AnalyzeRegionsTest("elt.3", "chr10a", 2, 10, "chr20c", 11, 19, false,
            new Fasta(fastaPath), new Fasta(fastaPath),
            RegionType.IDENTICAL),
        new AnalyzeRegionsTest("elt.4", "chr10c", 5, 13, "chr20d", 21, 29, false,
            new Fasta(fastaPath), new Fasta(fastaPath),
            RegionType.MISMATCHED_BASES),
        new AnalyzeRegionsTest("elt.5", "chr10d", 10, 18, "chr20e", 31, 39, true,
            new Fasta(fastaPath), new Fasta(fastaPath),
            RegionType.MISMATCHED_BASES),
        new AnalyzeRegionsTest("elt.6", "chr10d", 10, 18, "chr20f", 31, 39, false,
            new Fasta(fastaPath), new Fasta(fastaPath),
            RegionType.MISMATCHED_BASES),
        new AnalyzeRegionsTest("elt.7", "chr10d", 10, 18, "chr20x", 31, 39, true,
            new Fasta(fastaPath), new Fasta(fastaPath),
            RegionType.UNKNOWN_REGION_TYPE),
        new AnalyzeRegionsTest("elt.8", "chr10a", 2, 10, "chr20g", 12, 20, true,
            new Fasta(fastaPath), new Fasta(fastaPath),
            RegionType.UNKNOWN_REGION_TYPE),
        new AnalyzeRegionsTest("elt.8", "chr10a", 2, 10, "chr20h", 12, 20, true,
            new Fasta(fastaPath), new Fasta(fastaPath),
            RegionType.UNKNOWN_REGION_TYPE),
        new AnalyzeRegionsTest("elt.8", "chr10e", 2, 10, "chr20h", 12, 20, true,
            new Fasta(fastaPath), new Fasta(fastaPath),
            RegionType.UNKNOWN_REGION_TYPE),
        new AnalyzeRegionsTest("elt.8", "chr10f", 2, 10, "chr20i", 12, 20, true,
            new Fasta(fastaPath), new Fasta(fastaPath),
            RegionType.UNKNOWN_REGION_TYPE)
      };

      for (int i = 0; i < in.length; i++) {
        checkRegion(in[i]);
      }
    }

    @Test
    public void testJoinRegions() {
      List<GenomeRange> inQuery = new ArrayList<>();
      List<GenomeRange> inTarget = new ArrayList<>();

      inQuery.add(new GenomeRange("chr1", 1, 10, "elt.6", true));
      inQuery.add(new GenomeRange("chr2", 11, 20, "elt.5", true));
      inQuery.add(new GenomeRange("chr3", 21, 30, "elt.4", true));
      inQuery.add(new GenomeRange("chr4", 31, 40, "elt.3", true));
      inQuery.add(new GenomeRange("chr5", 41, 50, "elt.2", true));
      inQuery.add(new GenomeRange("chr6", 51, 60, "elt.1", true));

      inTarget.add(new GenomeRange("chr1", 1, 10, "elt.5", true));
      inTarget.add(new GenomeRange("chr2", 11, 20, "elt.3", true));
      inTarget.add(new GenomeRange("chr3", 21, 30, "elt.2", true));
      inTarget.add(new GenomeRange("chr4", 31, 40, "elt.4", true));
      inTarget.add(new GenomeRange("chr5", 41, 50, "elt.6", true));
      inTarget.add(new GenomeRange("chr6", 51, 60, "elt.1", true));

      List<HomologousRange> want = new ArrayList<>();
      want.add(
          newHRange(
            newRange("chr1", 1, 10),
            newRange("chr5", 41, 50), true));
      want.add(
          newHRange(
            newRange("chr2", 11, 20),
            newRange("chr1", 1, 10), true));
      want.add(
          newHRange(
            newRange("chr3", 21, 30),
            newRange("chr4", 31, 40), true));
      want.add(
          newHRange(
            newRange("chr4", 31, 40),
            newRange("chr2", 11, 20), true));
      want.add(
          newHRange(
            newRange("chr5", 41, 50),
            newRange("chr3", 21, 30), true));
      want.add(
          newHRange(
            newRange("chr6", 51, 60),
            newRange("chr6", 51, 60), true));

      List<HomologousRange> got = GenomeWarpSerial.joinRegions(inQuery, inTarget);
      assertTrue(GenomeWarpTestUtils.equivalentRangeWithoutType(got, want));
    }
  }

  /**
   * Tests for main GenomeWarp operation
   */
  private static final String QUERY_GENOME_PATH = "query.fasta";

  private static final String TARGET_GENOME_PATH = "target.fasta";

  private static final List<String> CALLSET_NAMES = new ArrayList<String>() {{
    add("myCallSetName1");
    add("myCallSetName2");
    add("myCallSetName3");
  }};

  private static final HomologousRange LARGEST_CHR1_HOMOLOG = newHRange(
      newRange("chr1", 1, 40), newRange("chr1_same", 11, 50), true, RegionType.IDENTICAL);

  private static final HomologousRange LARGEST_CHR2_REVCOMP_HOMOLOG = newHRange(
      newRange("chr2", 1, 43), newRange("chr2_revcomp", 10, 52), false, RegionType.IDENTICAL);

  private static final HomologousRange CHR2_MISMATCHED_BASES_HOMOLOG = newHRange(
      newRange("chr2", 0, 75), newRange("chr2_mismatched_bases", 0, 75), true,
      RegionType.MISMATCHED_BASES);

  private static final HomologousRange CHR2_CTG_DELETION_ALIGNMENT_HOMOLOG = newHRange(
      newRange("chr2", 1, 74), newRange("chr2_CTG_deletion", 11, 81), true,
      RegionType.ALIGNMENT_REQUIRED);

  private static final HomologousRange LARGEST_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG = newHRange(
      newRange("chr2", 1, 19), newRange("chr2_CTG_deletion", 11, 29), true,
      RegionType.IDENTICAL);

  private static final List<Variant> EMPTY_VARIANT_LIST = ImmutableList.of();

  private static final List<VariantCall> HOMOZYGOUS_REFERENCE =
      ImmutableList.of(
          VariantCall.newBuilder()
              .setCallSetName(CALLSET_NAMES.get(0))
              .addGenotype(0)
              .addGenotype(0)
              .build(),
          VariantCall.newBuilder()
              .setCallSetName(CALLSET_NAMES.get(1))
              .addGenotype(0)
              .addGenotype(0)
              .build(),
          VariantCall.newBuilder()
              .setCallSetName(CALLSET_NAMES.get(2))
              .addGenotype(0)
              .addGenotype(0)
              .build());
  private static final List<VariantCall> HETEROZYGOUS =
      ImmutableList.of(
          VariantCall.newBuilder()
              .setCallSetName(CALLSET_NAMES.get(0))
              .addGenotype(0)
              .addGenotype(1)
              .build());
  private static final List<VariantCall> HET2 =
      ImmutableList.of(
          VariantCall.newBuilder()
              .setCallSetName(CALLSET_NAMES.get(0))
              .addGenotype(1)
              .addGenotype(0)
              .build());
  private static final List<VariantCall> HOMOZYGOUS_ALTERNATE =
      ImmutableList.of(
          VariantCall.newBuilder()
              .setCallSetName(CALLSET_NAMES.get(0))
              .addGenotype(1)
              .addGenotype(1)
              .build());

  private static Variant createVariant(String chr, long start, String ref, String alt,
      List<VariantCall> calls) {
    return Variant.newBuilder().addFilter("PASS").setReferenceName(chr).setStart(start)
        .setEnd(start + ref.length()).setReferenceBases(ref).addAlternateBases(alt)
        .addAllCalls(calls).build();
  }

  @RunWith(Parameterized.class)
  public static final class PerformVariantCoordinateTransformTest {

    @Parameters
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            // TEST 0: Empty identical range should return the range with no variants in it.
            {
              LARGEST_CHR1_HOMOLOG,
              EMPTY_VARIANT_LIST,
              EMPTY_VARIANT_LIST
            },
            // TEST 1: Test simple transformations with different variant calls.
            {
              // Homozygous reference is not reported.
              LARGEST_CHR1_HOMOLOG,
              Arrays.asList(createVariant("chr1", 3, "G", "T", HOMOZYGOUS_REFERENCE)),
              Arrays.asList(createVariant("chr1_same", 13, "G", "T", HOMOZYGOUS_REFERENCE))
            },
            // TEST 2
            {
              LARGEST_CHR1_HOMOLOG,
              Arrays.asList(createVariant("chr1", 4, "T", "A", HETEROZYGOUS)),
              Arrays.asList(createVariant("chr1_same", 14, "T", "A", HETEROZYGOUS))
            },
            // TEST 3
            {
              LARGEST_CHR1_HOMOLOG,
              Arrays.<Variant>asList(createVariant("chr1", 5, "G", "GA", HOMOZYGOUS_ALTERNATE)),
              Arrays.<Variant>asList(createVariant("chr1_same", 15, "G", "GA",
                  HOMOZYGOUS_ALTERNATE))
            },
            // TEST 4
            {
              LARGEST_CHR1_HOMOLOG,
              Arrays.asList(
                  createVariant("chr1", 3, "G", "T", HOMOZYGOUS_REFERENCE),
                  createVariant("chr1", 5, "G", "GA", HOMOZYGOUS_ALTERNATE)),
              Arrays.asList(
                          createVariant("chr1_same", 13, "G", "T", HOMOZYGOUS_REFERENCE),
                          createVariant("chr1_same", 15, "G", "GA", HOMOZYGOUS_ALTERNATE))
            },
            // TEST 5
            {
              LARGEST_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              Arrays.asList(createVariant("chr2", 1, "ACT", "A", HETEROZYGOUS)),
              Arrays.asList(createVariant("chr2_CTG_deletion", 11, "ACT", "A", HETEROZYGOUS))
            },
            // TEST 6
            {
              LARGEST_CHR2_CTG_DELETION_IDENTICAL_HOMOLOG,
              Arrays.asList(createVariant("chr2", 1, "ACTG", "A", HETEROZYGOUS)),
              Arrays.asList(createVariant("chr2_CTG_deletion", 11, "A", "ACTG", HET2))
            },
            // TEST 7
            {
              LARGEST_CHR2_REVCOMP_HOMOLOG,
              EMPTY_VARIANT_LIST,
              EMPTY_VARIANT_LIST
            },
            // TEST 8
            {
              LARGEST_CHR2_REVCOMP_HOMOLOG,
              Arrays.asList(createVariant("chr2", 3, "T", "C", HOMOZYGOUS_REFERENCE)),
              Arrays.asList(createVariant("chr2_revcomp", 49, "A", "G", HOMOZYGOUS_REFERENCE))
            },
            // TEST 9
            {
              LARGEST_CHR2_REVCOMP_HOMOLOG,
              Arrays.asList(createVariant("chr2", 3, "T", "C", HETEROZYGOUS)),
              Arrays.asList(createVariant("chr2_revcomp", 49, "A", "G", HETEROZYGOUS))
            },
            // TEST 10: Empty alignment-required range -- not currently supported, should
            // be screened out.
            {
              CHR2_CTG_DELETION_ALIGNMENT_HOMOLOG,
              EMPTY_VARIANT_LIST,
              null
            },
            // TEST 11
            {
              CHR2_CTG_DELETION_ALIGNMENT_HOMOLOG,
              Arrays.asList(createVariant("chr2", 4, "G", "CA", HETEROZYGOUS)),
              null
            },
            // TEST 12: Mismatched bases in reference assembly should produce variants.
            {
              CHR2_MISMATCHED_BASES_HOMOLOG,
              EMPTY_VARIANT_LIST,
              Arrays.asList(
                  createVariant("chr2_mismatched_bases", 6, "A", "T", HOMOZYGOUS_ALTERNATE),
                  createVariant(
                      "chr2_mismatched_bases", 31, "T", "A", HOMOZYGOUS_ALTERNATE),
                  createVariant(
                      "chr2_mismatched_bases", 41, "T", "A", HOMOZYGOUS_ALTERNATE),
                  createVariant(
                      "chr2_mismatched_bases", 54, "T", "A", HOMOZYGOUS_ALTERNATE),
                  createVariant(
                      "chr2_mismatched_bases", 64, "T", "A", HOMOZYGOUS_ALTERNATE),
                  createVariant(
                      "chr2_mismatched_bases", 73, "T", "A", HOMOZYGOUS_ALTERNATE))
            },
            // TEST 13: Mismatched bases in reference assembly should produce variants
            {
              CHR2_MISMATCHED_BASES_HOMOLOG,
              EMPTY_VARIANT_LIST,
              Arrays.asList(
                  createVariant("chr2_mismatched_bases", 6, "A", "T", HOMOZYGOUS_ALTERNATE),
                  createVariant(
                      "chr2_mismatched_bases", 31, "T", "A", HOMOZYGOUS_ALTERNATE),
                  createVariant(
                      "chr2_mismatched_bases", 41, "T", "A", HOMOZYGOUS_ALTERNATE),
                  createVariant(
                      "chr2_mismatched_bases", 54, "T", "A", HOMOZYGOUS_ALTERNATE),
                  createVariant(
                      "chr2_mismatched_bases", 64, "T", "A", HOMOZYGOUS_ALTERNATE),
                  createVariant(
                      "chr2_mismatched_bases", 73, "T", "A", HOMOZYGOUS_ALTERNATE))
            },
          });
    }

    private HomologousRange inputRegion;

    private List<Variant> queryVariants;

    private List<Variant> pExpected;

    public PerformVariantCoordinateTransformTest(HomologousRange input1, List<Variant> input2,
        List<Variant> expected) {
      inputRegion = input1;
      queryVariants = input2;
      pExpected = expected;
    }

    @Test
    public void performVariantCoordinateTransform() throws Exception {
      Fasta queryFasta = new Fasta(GenomeWarpSerialTest.class.getClassLoader()
          .getResource(QUERY_GENOME_PATH).getFile());
      Fasta targetFasta = new Fasta(GenomeWarpSerialTest.class.getClassLoader()
          .getResource(TARGET_GENOME_PATH).getFile());

      List<Variant> actual = TransformVariantsInRange.transformQueryToTargetVariants(inputRegion,
          queryVariants, CALLSET_NAMES, queryFasta, targetFasta);

      assertTrue((actual == null && pExpected == null) || actual.equals(pExpected));
    }
  }
}
