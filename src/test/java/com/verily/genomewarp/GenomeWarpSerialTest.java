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
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.genomics.v1.Range;
import com.google.genomics.v1.Variant;
import com.google.genomics.v1.VariantCall;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.RegionType;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange.TargetStrand;
import com.verily.genomewarp.utils.Fasta;
import com.verily.genomewarp.utils.GenomeRange;
import com.verily.genomewarp.utils.GenomeRangeUtils;
import com.verily.genomewarp.utils.GenomeWarpTestUtils;
import htsjdk.variant.vcf.VCFFileReader;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
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
    private static final String VCF_PATH = "test.vcf";


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
      SortedMap<String, List<GenomeRange>> want = new TreeMap<String, List<GenomeRange>>() {{
        put("chr1", new ArrayList<GenomeRange>() {{
            add(new GenomeRange("chr1", 10, 24));
        }});
        put("chr2", new ArrayList<GenomeRange>() {{
          add(new GenomeRange("chr2", 5, 6));
        }});
        put("chr5", new ArrayList<GenomeRange>() {{
          add(new GenomeRange("chr5", 15, 20));
        }});
        put("chr6", new ArrayList<GenomeRange>() {{
          add(new GenomeRange("chr6", 30, 35));
        }});
        put("chr7", new ArrayList<GenomeRange>() {{
          add(new GenomeRange("chr7", 40, 44));
          add(new GenomeRange("chr7", 45, 50));
        }});
        put("chr8", new ArrayList<GenomeRange>() {{
          add(new GenomeRange("chr8", 56, 60));
          add(new GenomeRange("chr8", 61, 66));
          add(new GenomeRange("chr8", 67, 72));
        }});
      }};

      // Get test result
      SortedMap<String, List<GenomeRange>> got =
          GenomeRangeUtils.splitAtNonDNA(testFasta, testBED);

      assertEquals(got.size(), want.size());
      for(String key: got.keySet()) {
        assertTrue(GenomeWarpTestUtils.equivalentRanges(got.get(key), want.get(key)));
      }
    }

    @Test
    public void testGenerateBEDFromVCF() {
      String file = GenomeWarpSerialTest.class.getClassLoader().getResource(VCF_PATH).getFile();
      VCFFileReader testVcf = new VCFFileReader(new File(file), false);

      SortedMap<String, List<GenomeRange>> want = new TreeMap<String, List<GenomeRange>>() {{
        put("chr1", new ArrayList<GenomeRange>() {{
          add(new GenomeRange("chr1", 49, 50));
          add(new GenomeRange("chr1", 51, 52));
          add(new GenomeRange("chr1", 119, 122));
          add(new GenomeRange("chr1", 136, 137));
          add(new GenomeRange("chr1", 189, 190));
        }});
        put("chr2", new ArrayList<GenomeRange>() {{
          add(new GenomeRange("chr2", 139, 141));
        }});
      }};

      SortedMap<String, List<GenomeRange>> got = GenomeRangeUtils.generateBEDFromVCF(testVcf);
      assertEquals(got.size(), want.size());
      for(String key: got.keySet()) {
        assertTrue(GenomeWarpTestUtils.equivalentRanges(got.get(key), want.get(key)));
      }
    }

    private void checkSplitRegions(GenomeRange in,
        GenomeRange[] want) {
      List<GenomeRange> wantList;
      wantList = Arrays.asList(want);

      List<GenomeRange> got =
          GenomeRangeUtils.splitRegion(in, 10);

      assertTrue(GenomeWarpTestUtils.equivalentRanges(got, wantList));
    }

    @Test
    public void testSplitRegions() {
      GenomeRange in1 = new GenomeRange("chr1", 4, 10);
      GenomeRange[] want1 = {
          new GenomeRange("chr1", 4, 10)
      };
      checkSplitRegions(in1, want1);

      GenomeRange in2 = new GenomeRange("chr1", 0, 10);
      GenomeRange[] want2 = {
          new GenomeRange("chr1", 0, 10)
      };
      checkSplitRegions(in2, want2);

      GenomeRange in3 = new GenomeRange("chr1", 18, 45);
      GenomeRange[] want3 = {
          new GenomeRange("chr1", 18, 28),
          new GenomeRange("chr1", 28, 38),
          new GenomeRange("chr1", 38, 45),
      };
      checkSplitRegions(in3, want3);
    }

    @Test
    public void testFilterOutNonCoveredVariants() {
      GenomeRange[] queryRegions = {
          new GenomeRange("chr1", 4, 10),
          new GenomeRange("chr1", 12345, 12545),
          new GenomeRange("chr1", 15123, 15130),
          new GenomeRange("chr1", 20000, 40000),
      };

      GenomeRange[] variantRegions = {
          new GenomeRange("chr1", 0, 4),
          new GenomeRange("chr1", 10, 15),
          new GenomeRange("chr1", 123, 345),
          new GenomeRange("chr1", 12340, 12348),
          new GenomeRange("chr1", 12345, 12350),
          new GenomeRange("chr1", 12400, 12407),
          new GenomeRange("chr1", 12538, 12540),
          new GenomeRange("chr1", 12538, 12545),
          new GenomeRange("chr1", 12538, 12560),
          new GenomeRange("chr1", 15120, 15140),
          new GenomeRange("chr1", 15123, 15130),
          new GenomeRange("chr1", 100000, 100005),
      };

      GenomeRange[] want = {
          new GenomeRange("chr1", 12345, 12350),
          new GenomeRange("chr1", 12400, 12407),
          new GenomeRange("chr1", 12538, 12540),
          new GenomeRange("chr1", 12538, 12545),
          new GenomeRange("chr1", 15123, 15130),
      };

      List<GenomeRange> got = GenomeRangeUtils.filterOutNotCoveredVariants(
          Arrays.asList(queryRegions), Arrays.asList(variantRegions));
      assertTrue(GenomeWarpTestUtils.equivalentRanges(Arrays.asList(want), got));
    }

    @Test
    public void testGeneratePaddedGenomeRanges() {
      GenomeRange[] in = {
          new GenomeRange("chr1", 12345, 12350),
          new GenomeRange("chr1", 12400, 12407),
          new GenomeRange("chr1", 12538, 12540),
          new GenomeRange("chr1", 12538, 12545),
          new GenomeRange("chr1", 15123, 15130),
      };

      GenomeRange[] want = {
          new GenomeRange("chr1", 12340, 12355),
          new GenomeRange("chr1", 12395, 12412),
          new GenomeRange("chr1", 12533, 12545),
          new GenomeRange("chr1", 12533, 12550),
          new GenomeRange("chr1", 15118, 15135),
      };

      List<GenomeRange> got = GenomeRangeUtils.generatePaddedGenomeRanges(Arrays.asList(in));
      assertTrue(GenomeWarpTestUtils.equivalentRanges(Arrays.asList(want), got));
    }

    private void checkMergeRegionsFromQueryBEDAndVariants(
        GenomeRange[] queryBED, GenomeRange[] vcfBED, GenomeRange[] want) {
      List<GenomeRange> wantList = Arrays.asList(want);

      List<GenomeRange> got = GenomeRangeUtils.mergeRegionsFromQueryBEDAndVariants(
              Arrays.asList(queryBED), Arrays.asList(vcfBED), 10);

      assertTrue(GenomeWarpTestUtils.equivalentRanges(got, wantList));
    }

    @Test
    public void testMergeRegionsFromQueryBEDAndVariants() {
      GenomeRange[] queryBED1 = {
        new GenomeRange("chr1", 90, 95),
        new GenomeRange("chr1", 100, 150),
      };
      GenomeRange[] vcfBED1 = {
        new GenomeRange("chr1", 108, 125),
      };
      GenomeRange[] want1 = {
        new GenomeRange("chr1", 90, 95),
        new GenomeRange("chr1", 100, 108),
        new GenomeRange("chr1", 108, 125),
        new GenomeRange("chr1", 125, 135),
        new GenomeRange("chr1", 135, 145),
        new GenomeRange("chr1", 145, 150),
      };

      checkMergeRegionsFromQueryBEDAndVariants(queryBED1, vcfBED1, want1);

      GenomeRange[] queryBED2 = {
        new GenomeRange("chr1", 90, 95),
        new GenomeRange("chr1", 100, 150),
      };
      GenomeRange[] vcfBED2 = {
      };
      GenomeRange[] want2 = {
        new GenomeRange("chr1", 90, 95),
        new GenomeRange("chr1", 100, 110),
        new GenomeRange("chr1", 110, 120),
        new GenomeRange("chr1", 120, 130),
        new GenomeRange("chr1", 130, 140),
        new GenomeRange("chr1", 140, 150),
      };

      checkMergeRegionsFromQueryBEDAndVariants(queryBED2, vcfBED2, want2);

      GenomeRange[] queryBED3 = {
      };
      GenomeRange[] vcfBED3 = {
        new GenomeRange("chr1", 108, 125),
      };
      GenomeRange[] want3 = {
      };

      checkMergeRegionsFromQueryBEDAndVariants(queryBED3, vcfBED3, want3);

      GenomeRange[] queryBED4 = {
        new GenomeRange("chr1", 90, 95),
        new GenomeRange("chr1", 120, 150),
        new GenomeRange("chr1", 152, 160),
        new GenomeRange("chr1", 170, 180),
        new GenomeRange("chr1", 2000, 2005),
      };
      GenomeRange[] vcfBED4 = {
        new GenomeRange("chr1", 88, 97), // includes confident region
        new GenomeRange("chr1", 115, 128), // begins before confident region
        new GenomeRange("chr1", 145, 155), // includes fragments from 2 confident region
        new GenomeRange("chr1", 170, 175), // simple
        new GenomeRange("chr1", 178, 188), // the same conf region as previos variant
      };
      GenomeRange[] want4 = {
        new GenomeRange("chr1", 90, 95),
        new GenomeRange("chr1", 120, 128),
        new GenomeRange("chr1", 128, 138),
        new GenomeRange("chr1", 138, 145),
        new GenomeRange("chr1", 145, 150),
        new GenomeRange("chr1", 152, 155),
        new GenomeRange("chr1", 155, 160),
        new GenomeRange("chr1", 170, 175),
        new GenomeRange("chr1", 175, 178),
        new GenomeRange("chr1", 178, 180),
        new GenomeRange("chr1", 2000, 2005),
      };

      checkMergeRegionsFromQueryBEDAndVariants(queryBED4, vcfBED4, want4);
    }

    @Test
    public void testMergeOverlaps() {
      GenomeRange[] in = {
        new GenomeRange("chr1", 90, 95),
              new GenomeRange("chr1", 100, 155),
              new GenomeRange("chr1", 150, 200),
              new GenomeRange("chr1", 325, 335),
              new GenomeRange("chr1", 335, 350),
              new GenomeRange("chr2", 10, 50),
              new GenomeRange("chr2", 80, 90),
              new GenomeRange("chr2", 85, 120),
              new GenomeRange("chr2", 120, 150),
      };
      GenomeRange[] want = {
              new GenomeRange("chr1", 90, 95),
              new GenomeRange("chr1", 100, 200),
              new GenomeRange("chr1", 325, 350),
              new GenomeRange("chr2", 10, 50),
              new GenomeRange("chr2", 80, 150),
      };

      List<GenomeRange> got = GenomeRangeUtils.mergeOverlaps(Arrays.asList(in));

      assertTrue(GenomeWarpTestUtils.equivalentRanges(got, Arrays.asList(want)));
    }

    @Test
    public void testMergeOverlapsUnsorted() {
      thrown.expect(RuntimeException.class);
      thrown.expectMessage("input regions are not sorted by position");
      GenomeRange[] in = {
              new GenomeRange("chr1", 108, 125),
              new GenomeRange("chr1", 108, 115),
              new GenomeRange("chr2", 10, 15),
      };
      GenomeRangeUtils.mergeOverlaps(Arrays.asList(in));
    }

    private void checkOmitOverlaps(GenomeRange[] in,
        GenomeRange[] want) {
      List<GenomeRange> inList, wantList;
      inList = Arrays.asList(in);
      wantList = Arrays.asList(want);

      List<GenomeRange> gotList =
          GenomeRangeUtils.omitOverlap(inList);

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
      };
      GenomeRange[] want3 = {
          new GenomeRange("chr1", 50, 56),
          new GenomeRange("chr1", 60, 66),
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
      };
      GenomeRange[] want4 = {
          new GenomeRange("chr1", 60, 66),
      };
      checkOmitOverlaps(in4, want4);
    }

    @Rule
    public ExpectedException thrown = ExpectedException.none();

    private void checkOmitOverlapsBad(GenomeRange[] in) {
      List<GenomeRange> inList = Arrays.asList(in);

      GenomeRangeUtils.omitOverlap(inList);
    }

    @Test
    public void testOmitOverlapUnsortedPos() {
      thrown.expect(RuntimeException.class);
      thrown.expectMessage("output BED of liftover is not sorted by position");
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
    public void testOmitOverlapDifferentChr() {
      thrown.expect(RuntimeException.class);
      thrown.expectMessage("found ranges from different chromosomes");
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

    @Test
    public void testSimplifiedRegionsPreprocessing() {
      SortedMap<String, List<GenomeRange>> in = new TreeMap<String, List<GenomeRange>>() {{
          put("chr1", new ArrayList<GenomeRange>() {{
              add(new GenomeRange("chr1", 30, 60));
          }});
          put("chr2", new ArrayList<GenomeRange>() {{
              add(new GenomeRange("chr2", 50, 60));
              add(new GenomeRange("chr2", 1000, 1020));
              add(new GenomeRange("chr2", 5, 6));
          }});
      }};
      GenomeRange[] want = {
          new GenomeRange("chr1", 30, 60),
          new GenomeRange("chr2", 5, 6),
          new GenomeRange("chr2", 50, 60),
          new GenomeRange("chr2", 1000, 1020),
      };

      List<GenomeRange> got = GenomeRangeUtils.generateQueryBEDWithSimplifiedPreprocessing(in);
      assertTrue(GenomeWarpTestUtils.equivalentRanges(got, Arrays.asList(want)));
      assertTrue(!got.get(0).getName().isEmpty());
    }

    @Test
    public void testImprovedRegionsPreprocessing() {
      SortedMap<String, List<GenomeRange>> in = new TreeMap<String, List<GenomeRange>>() {{
        put("chr1", new ArrayList<GenomeRange>() {{
          add(new GenomeRange("chr1", 30, 60));
        }});
        put("chr2", new ArrayList<GenomeRange>() {{
          add(new GenomeRange("chr2", 50, 60));
          add(new GenomeRange("chr2", 1000, 1020));
          add(new GenomeRange("chr2", 5, 6));
        }});
      }};

      // Variants: chr1 - 49:50, 51:52, 119:122, 136:137, 189:190
      // chr2 - 139:141

      GenomeRange[] want = {
              new GenomeRange("chr1", 30, 40),
              new GenomeRange("chr1", 40, 44),
              new GenomeRange("chr1", 44, 57),
              new GenomeRange("chr1", 57, 60),
              new GenomeRange("chr2", 5, 6),
              new GenomeRange("chr2", 50, 60),
              new GenomeRange("chr2", 1000, 1010),
              new GenomeRange("chr2", 1010, 1020),
      };

      List<GenomeRange> got = GenomeRangeUtils.generateQueryBEDWithImprovedPreprocessing(
              in, GenomeWarpSerialTest.class.getClassLoader().getResource(VCF_PATH).getFile(), 10);
      assertTrue(GenomeWarpTestUtils.equivalentRanges(got, Arrays.asList(want)));
      assertTrue(!got.get(0).getName().isEmpty());
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

      List<HomologousRange> got = GenomeRangeUtils.joinRegions(inQuery, inTarget);
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
