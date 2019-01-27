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

import static com.google.common.base.Preconditions.checkArgument;

import com.google.genomics.v1.Variant;
import com.google.genomics.v1.VariantCall;
import com.google.protobuf.ListValue;
import com.google.protobuf.Value;
import com.verily.genomewarp.HomologousRangeOuterClass.HomologousRange;
import htsjdk.variant.vcf.VCFConstants;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Utility functions and variables for various GenomeWarp tests
 */
public class GenomeWarpTestUtils {

  private static Value valueFromObject(Object v) {
    if (v instanceof Boolean) {
      return Value.newBuilder().setBoolValue((Boolean) v).build();
    }
    if (v instanceof Integer) {
      return Value.newBuilder().setNumberValue((Integer) v).build();
    }
    if (v instanceof Double) {
      return Value.newBuilder().setNumberValue((Double) v).build();
    }

    return Value.newBuilder().setStringValue((String) v).build();
  }

  private static ListValue listFromObject(Object obj) {
    ListValue.Builder lvBuilder = ListValue.newBuilder();
    if (!(obj instanceof List)) {
      lvBuilder.addValues(valueFromObject(obj));
    } else {
      List<Object> objList = (List<Object>) obj;

      for (int i = 0; i < objList.size(); i++) {
        lvBuilder.addValues(valueFromObject(objList.get(i)));
      }
    }

    return lvBuilder.build();
  }

  public static Variant makeVariant(String chr, long start, String ref, String... alts) {
    return Variant.newBuilder()
        .addFilter(VCFConstants.PASSES_FILTERS_v4)
        .setReferenceName(chr)
        .setStart(start)
        .setEnd(start + ref.length())
        .setReferenceBases(ref)
        .addAllAlternateBases(Arrays.asList(alts))
        .build();
  }

  public static Variant makeVariant(
      String chr, long start, String ref, String alt, List<VariantCall> calls) {
    return Variant.newBuilder()
        .addFilter(VCFConstants.PASSES_FILTERS_v4)
        .setReferenceName(chr)
        .setStart(start)
        .setEnd(start + ref.length())
        .setReferenceBases(ref)
        .addAlternateBases(alt)
        .addAllCalls(calls)
        .build();
  }

  private static VariantCall makeVariantCall(
      String name, int[] g, double[] gl, boolean isPhased, Map<String, ListValue> info) {
    VariantCall.Builder vcBuilder = VariantCall.newBuilder().setCallSetName(name);
    if (g != null) {
      for (int i = 0; i < g.length; i++) {
        vcBuilder.addGenotype(g[i]);
      }
      if (isPhased && g.length > 1) {
        vcBuilder.setPhaseset("*");
      }
    }
    if (gl != null) {
      for (int i = 0; i < gl.length; i++) {
        vcBuilder.addGenotypeLikelihood(gl[i]);
      }
    }
    if (info != null) {
      vcBuilder.getMutableInfo().putAll(info);
    }

    return vcBuilder.build();
  }

  public static boolean equivalentRangeWithoutType(List<HomologousRange> listOne,
      List<HomologousRange> listTwo) {
    if (listOne.size() != listTwo.size()) {
      return false;
    }
    Collections.sort(listOne, new GenomeWarpUtils.HomologousRangeComparator());
    Collections.sort(listTwo, new GenomeWarpUtils.HomologousRangeComparator());
    int index = 0;
    for (HomologousRange ruleOne : listOne) {
      HomologousRange ruleTwo = listTwo.get(index++);
      if (!ruleOne.getQueryRange().equals(ruleTwo.getQueryRange())) {
        return false;
      }
      if (!ruleOne.getTargetRange().equals(ruleTwo.getTargetRange())) {
        return false;
      }
      if (!ruleOne.getTargetStrand().equals(ruleTwo.getTargetStrand())) {
        return false;
      }
    }

    return true;
  }

  /**
   * Returns true if and only if the items in each list are identical.
   *
   * This function does not require the ordering of the ranges to be the same,
   * but does require the counts of each item to be the same.
   *
   * @param listOne a list of GenomeRanges to compare.
   * @param listTwo another list of GenomeRanges to compare.
   * @returns true iff the lists contain the same elements.
   */
  public static boolean equivalentRanges(List<GenomeRange> listOne,
      List<GenomeRange> listTwo) {
    if (listOne.size() != listTwo.size()) {
      return false;
    }
    Collections.sort(listOne);
    Collections.sort(listTwo);
    int index = 0;
    for (GenomeRange range : listOne) {
      if (!range.isSameInterval(listTwo.get(index++))) {
        return false;
      }
    }

    return true;
  }

  /**
   * Returns true if and only if the items in each Map are identical.
   *
   * This requires the same chromosomes to be present in each Map, and the per-chromosome
   * ranges to contain the same items. It does not require the ordering of the ranges to be
   * the same within chromosomes.
   *
   * @param first a Map from chromosome -> GenomeRange list to compare.
   * @param second another Map from chromosome -> GenomeRange list to compare.
   * @returns true iff the maps contain the same elements.
   */
  public static boolean equivalentGenomeRanges(Map<String, List<GenomeRange>> first,
      Map<String, List<GenomeRange>> second) {
    if (first.size() != second.size()) {
      return false;
    }
    for (String key : first.keySet()) {
      if (!second.containsKey(key)) {
        return false;
      }
      if (!equivalentRanges(first.get(key), second.get(key))) {
        return false;
      }
    }
    return true;
  }

  /*
   * Build truth variants
   *
   * Note that currently, this tool only supports up to VCF v4.1. As a result,
   * variants 6, 9, 10, and 11 are not added to our test set and are omitted from
   * valid-4.1.vcf.
   */
  public static List<Variant> getTruth(boolean supportsV42) {
    List<Variant> truth = new ArrayList<>();

    // Variant 1
    Variant.Builder variantBuilder1 =
        Variant.newBuilder()
            .setReferenceName("19")
            .setStart(14369)
            .setEnd(14370)
            .addNames("rs6054257")
            .setReferenceBases("G")
            .addAlternateBases("A")
            .setQuality(29)
            .addFilter(VCFConstants.PASSES_FILTERS_v4);
    variantBuilder1
        .getMutableInfo()
        .putAll(
            new HashMap<String, ListValue>() {
              {
                put("NS", listFromObject(3));
                put("DP", listFromObject(14));
                put("AF", listFromObject(0.5));
                put("DB", listFromObject(true));
                put("H2", listFromObject(true));
              }
            });
    variantBuilder1.addCalls(
        makeVariantCall(
            "NA00001",
            new int[] {0, 0},
            null,
            true,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(99));
                put("DP", listFromObject(1));
                put("HQ", listFromObject(Arrays.asList(51, 51)));
              }
            }));
    variantBuilder1.addCalls(
        makeVariantCall(
            "NA00002",
            new int[] {1, 0},
            null,
            true,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(48));
                put("DP", listFromObject(8));
                put("HQ", listFromObject(Arrays.asList(51, 51)));
              }
            }));
    variantBuilder1.addCalls(
        makeVariantCall(
            "NA00003",
            new int[] {1, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(43));
                put("DP", listFromObject(5));
              }
            }));
    truth.add(variantBuilder1.build());

    // Variant 2
    Variant.Builder variantBuilder2 =
        Variant.newBuilder()
            .setReferenceName("20")
            .setStart(17329)
            .setEnd(17330)
            .setReferenceBases("T")
            .addAlternateBases("A")
            .setQuality(3)
            .addFilter("q10");
    variantBuilder2
        .getMutableInfo()
        .putAll(
            new HashMap<String, ListValue>() {
              {
                put("NS", listFromObject(3));
                put("DP", listFromObject(11));
                put("AF", listFromObject(0.017));
              }
            });
    variantBuilder2.addCalls(
        makeVariantCall(
            "NA00001",
            new int[] {0, 0},
            null,
            true,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(49));
                put("DP", listFromObject(3));
                put("HQ", listFromObject(Arrays.asList(58, 50)));
                put("FT", listFromObject("q10"));
              }
            }));
    variantBuilder2.addCalls(
        makeVariantCall(
            "NA00002",
            new int[] {0, 1},
            null,
            true,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(3));
                put("DP", listFromObject(5));
                put("HQ", listFromObject(Arrays.asList((Object) VCFConstants.MISSING_VALUE_v4,
                    (Object) Integer.valueOf(3))));
                put("FT", listFromObject("s50"));
              }
            }));
    variantBuilder2.addCalls(
        makeVariantCall(
            "NA00003",
            new int[] {0, 0},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(41));
                put("DP", listFromObject(3));
                put("FT", listFromObject("q10"));
              }
            }));
    truth.add(variantBuilder2.build());

    // Variant 3
    Variant.Builder variantBuilder3 =
        Variant.newBuilder()
            .setReferenceName("20")
            .setStart(1110695)
            .setEnd(1110696)
            .addNames("rs6040355")
            .setReferenceBases("A")
            .addAlternateBases("G")
            .addAlternateBases("T")
            .setQuality(67)
            .addFilter(VCFConstants.PASSES_FILTERS_v4);
    variantBuilder3
        .getMutableInfo()
        .putAll(
            new HashMap<String, ListValue>() {
              {
                put("NS", listFromObject(2));
                put("DP", listFromObject(10));
                put("AF", listFromObject(Arrays.asList(0.333, 0.667)));
                put("AA", listFromObject("T"));
                put("DB", listFromObject(true));
              }
            });
    variantBuilder3.addCalls(
        makeVariantCall(
            "NA00001",
            new int[] {1, 2},
            null,
            true,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(21));
                put("DP", listFromObject(6));
                put("HQ", listFromObject(Arrays.asList(23, 27)));
              }
            }));
    variantBuilder3.addCalls(
        makeVariantCall(
            "NA00002",
            new int[] {2, 1},
            null,
            true,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(2));
                put("DP", listFromObject(0));
                put("HQ", listFromObject(Arrays.asList(18, 2)));
              }
            }));
    variantBuilder3.addCalls(
        makeVariantCall(
            "NA00003",
            new int[] {2, 2},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(35));
                put("DP", listFromObject(4));
              }
            }));
    truth.add(variantBuilder3.build());

    // Variant 4
    Variant.Builder variantBuilder4 =
        Variant.newBuilder()
            .setReferenceName("20")
            .setStart(1230236)
            .setEnd(1230237)
            .setReferenceBases("T")
            .setQuality(47)
            .addFilter(VCFConstants.PASSES_FILTERS_v4);
    variantBuilder4
        .getMutableInfo()
        .putAll(
            new HashMap<String, ListValue>() {
              {
                put("NS", listFromObject(3));
                put("DP", listFromObject(13));
                put("AA", listFromObject("T"));
              }
            });
    variantBuilder4.addCalls(
        makeVariantCall(
            "NA00001",
            new int[] {0, 0},
            null,
            true,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(54));
                put("DP", listFromObject(7));
                put("HQ", listFromObject(Arrays.asList(56, 60)));
              }
            }));
    variantBuilder4.addCalls(
        makeVariantCall(
            "NA00002",
            new int[] {0, 0},
            null,
            true,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(48));
                put("DP", listFromObject(4));
                put("HQ", listFromObject(Arrays.asList(51, 51)));
              }
            }));
    variantBuilder4.addCalls(
        makeVariantCall(
            "NA00003",
            new int[] {0, 0},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(61));
                put("DP", listFromObject(2));
              }
            }));
    truth.add(variantBuilder4.build());

    // Variant 5
    Variant.Builder variantBuilder5 =
        Variant.newBuilder()
            .setReferenceName("20")
            .setStart(1234566)
            .setEnd(1234569)
            .addNames("microsat1")
            .setReferenceBases("GTC")
            .addAlternateBases("G")
            .addAlternateBases("GTCTC")
            .setQuality(50)
            .addFilter(VCFConstants.PASSES_FILTERS_v4);
    variantBuilder5
        .getMutableInfo()
        .putAll(
            new HashMap<String, ListValue>() {
              {
                put("NS", listFromObject(3));
                put("DP", listFromObject(9));
                put("AA", listFromObject("G"));
              }
            });
    variantBuilder5.addCalls(
        makeVariantCall(
            "NA00001",
            new int[] {0, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(35));
                put("DP", listFromObject(4));
              }
            }));
    variantBuilder5.addCalls(
        makeVariantCall(
            "NA00002",
            new int[] {0, 2},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(17));
                put("DP", listFromObject(2));
              }
            }));
    variantBuilder5.addCalls(
        makeVariantCall(
            "NA00003",
            new int[] {1, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(40));
                put("DP", listFromObject(3));
              }
            }));
    truth.add(variantBuilder5.build());

    // Variant 6 -- REQUIRES VCF 4.2
    // 20	2234567	.	C	[13:123457[ACGC	50	PASS	SVTYPE=BND;NS=3;DP=9;AA=G	GT:GQ:DP	0/1:35:4	0/1:17:2	1/1:40:3
    Variant.Builder variantBuilder6 =
        Variant.newBuilder()
            .setReferenceName("20")
            .setStart(2234566)
            .setEnd(2234567)
            .setReferenceBases("C")
            .addAlternateBases("[13:123457[ACGC")
            .setQuality(50)
            .addFilter(VCFConstants.PASSES_FILTERS_v4);
    variantBuilder6
        .getMutableInfo()
        .putAll(
            new HashMap<String, ListValue>() {
              {
                put("SVTYPE", listFromObject("BND"));
                put("NS", listFromObject(3));
                put("DP", listFromObject(9));
                put("AA", listFromObject("G"));
              }
            });
    variantBuilder6.addCalls(
        makeVariantCall(
            "NA00001",
            new int[] {0, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(35));
                put("DP", listFromObject(4));
              }
            }));
    variantBuilder6.addCalls(
        makeVariantCall(
            "NA00002",
            new int[] {0, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(17));
                put("DP", listFromObject(2));
              }
            }));
    variantBuilder6.addCalls(
        makeVariantCall(
            "NA00003",
            new int[] {1, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(40));
                put("DP", listFromObject(3));
              }
            }));
    if (supportsV42) {
      truth.add(variantBuilder6.build());
    }

    // Variant 7
    Variant.Builder variantBuilder7 =
        Variant.newBuilder()
            .setReferenceName("20")
            .setStart(2234567)
            .setEnd(2234568)
            .setReferenceBases("C")
            .addAlternateBases("TC")
            .setQuality(50)
            .addFilter(VCFConstants.PASSES_FILTERS_v4);
    variantBuilder7
        .getMutableInfo()
        .putAll(
            new HashMap<String, ListValue>() {
              {
                put("SVTYPE", listFromObject("BND"));
                put("NS", listFromObject(3));
                put("DP", listFromObject(9));
                put("AA", listFromObject("G"));
              }
            });
    variantBuilder7.addCalls(
        makeVariantCall(
            "NA00001",
            new int[] {0, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(35));
                put("DP", listFromObject(4));
              }
            }));
    variantBuilder7.addCalls(
        makeVariantCall(
            "NA00002",
            new int[] {0, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(17));
                put("DP", listFromObject(2));
              }
            }));
    variantBuilder7.addCalls(
        makeVariantCall(
            "NA00003",
            new int[] {1, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(40));
                put("DP", listFromObject(3));
              }
            }));
    truth.add(variantBuilder7.build());

    // Variant 8
    Variant.Builder variantBuilder8 =
        Variant.newBuilder()
            .setReferenceName("20")
            .setStart(2234568)
            .setEnd(2234569)
            .setReferenceBases("C")
            .addAlternateBases("CT")
            .setQuality(50)
            .addFilter(VCFConstants.PASSES_FILTERS_v4);
    variantBuilder8
        .getMutableInfo()
        .putAll(
            new HashMap<String, ListValue>() {
              {
                put("SVTYPE", listFromObject("BND"));
                put("NS", listFromObject(3));
                put("DP", listFromObject(9));
                put("AA", listFromObject("G"));
              }
            });
    variantBuilder8.addCalls(
        makeVariantCall(
            "NA00001",
            new int[] {0, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(35));
                put("DP", listFromObject(4));
              }
            }));
    variantBuilder8.addCalls(
        makeVariantCall(
            "NA00002",
            new int[] {0, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(17));
                put("DP", listFromObject(2));
              }
            }));
    variantBuilder8.addCalls(
        makeVariantCall(
            "NA00003",
            new int[] {-1, -1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(40));
                put("DP", listFromObject(3));
              }
            }));
    truth.add(variantBuilder8.build());

    // Variant 9 -- REQUIRES VCF 4.2
    // 20	3234569	.	C	<SYMBOLIC>	50	PASS	END=3235677;NS=3;DP=9;AA=G	GT:GQ:DP	0/1:35:4	0/1:17:2	1/1:40:3
    Variant.Builder variantBuilder9 =
        Variant.newBuilder()
            .setReferenceName("20")
            .setStart(2234568)
            .setEnd(2234569)
            .setReferenceBases("C")
            .addAlternateBases("<SYMBOLIC>")
            .setQuality(50)
            .addFilter(VCFConstants.PASSES_FILTERS_v4);
    variantBuilder9
        .getMutableInfo()
        .putAll(
            new HashMap<String, ListValue>() {
              {
                put("END", listFromObject(3235677));
                put("NS", listFromObject(3));
                put("DP", listFromObject(9));
                put("AA", listFromObject("G"));
              }
            });
    variantBuilder9.addCalls(
        makeVariantCall(
            "NA00001",
            new int[] {0, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(35));
                put("DP", listFromObject(4));
              }
            }));
    variantBuilder9.addCalls(
        makeVariantCall(
            "NA00002",
            new int[] {0, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(17));
                put("DP", listFromObject(2));
              }
            }));
    variantBuilder9.addCalls(
        makeVariantCall(
            "NA00003",
            new int[] {1, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(40));
                put("DP", listFromObject(3));
              }
            }));
    if (supportsV42) {
      truth.add(variantBuilder9.build());
    }

    // Variant 10 -- REQUIRES VCF 4.2
    // 20	4234569	.	N	.[13:123457[	50	PASS	SVTYPE=BND;NS=3;DP=9;AA=G	GT:GQ:DP	0/1:35:4	0/1:17:2	./.:40:3
    Variant.Builder variantBuilder10 =
        Variant.newBuilder()
            .setReferenceName("20")
            .setStart(2234568)
            .setEnd(2234569)
            .setReferenceBases("N")
            .addAlternateBases(".[13:123457[")
            .setQuality(50)
            .addFilter(VCFConstants.PASSES_FILTERS_v4);
    variantBuilder10
        .getMutableInfo()
        .putAll(
            new HashMap<String, ListValue>() {
              {
                put("SVTYPE", listFromObject("BND"));
                put("NS", listFromObject(3));
                put("DP", listFromObject(9));
                put("AA", listFromObject("G"));
              }
            });
    variantBuilder10.addCalls(
        makeVariantCall(
            "NA00001",
            new int[] {0, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(35));
                put("DP", listFromObject(4));
              }
            }));
    variantBuilder10.addCalls(
        makeVariantCall(
            "NA00002",
            new int[] {0, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(17));
                put("DP", listFromObject(2));
              }
            }));
    variantBuilder10.addCalls(
        makeVariantCall(
            "NA00003",
            new int[] {-1, -1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(40));
                put("DP", listFromObject(3));
              }
            }));
    if (supportsV42) {
      truth.add(variantBuilder10.build());
    }

    // Variant 11 -- REQUIRES VCF 4.2
    // 20	5234569	.	N	[13:123457[.	50	PASS	SVTYPE=BND;NS=3;DP=9;AA=G	GT:GQ:DP	0/1:35:4	0/1:17:2	1/1:40:3
    Variant.Builder variantBuilder11 =
        Variant.newBuilder()
            .setReferenceName("20")
            .setStart(5234568)
            .setEnd(5234569)
            .setReferenceBases("N")
            .addAlternateBases("[13:123457[.")
            .setQuality(50)
            .addFilter(VCFConstants.PASSES_FILTERS_v4);
    variantBuilder11
        .getMutableInfo()
        .putAll(
            new HashMap<String, ListValue>() {
              {
                put("SVTYPE", listFromObject("BND"));
                put("NS", listFromObject(3));
                put("DP", listFromObject(9));
                put("AA", listFromObject("G"));
              }
            });
    variantBuilder11.addCalls(
        makeVariantCall(
            "NA00001",
            new int[] {0, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(35));
                put("DP", listFromObject(4));
              }
            }));
    variantBuilder11.addCalls(
        makeVariantCall(
            "NA00002",
            new int[] {0, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(17));
                put("DP", listFromObject(2));
              }
            }));
    variantBuilder11.addCalls(
        makeVariantCall(
            "NA00003",
            new int[] {1, 1},
            null,
            false,
            new HashMap<String, ListValue>() {
              {
                put("GQ", listFromObject(40));
                put("DP", listFromObject(3));
              }
            }));
    if (supportsV42) {
      truth.add(variantBuilder11.build());
    }

    // Variant 12
    Variant.Builder variantBuilder12 =
        Variant.newBuilder()
            .setReferenceName("Y")
            .setStart(17329)
            .setEnd(17330)
            .setReferenceBases("T")
            .addAlternateBases("A")
            .setQuality(3)
            .addFilter("q10");
    variantBuilder12
        .getMutableInfo()
        .putAll(
            new HashMap<String, ListValue>() {
              {
                put("NS", listFromObject(3));
                put("DP", listFromObject(11));
                put("AF", listFromObject(0.017));
              }
            });
    variantBuilder12.addCalls(
        makeVariantCall("NA00001", new int[] {0}, new double[] {-49, 0}, false, null));
    variantBuilder12.addCalls(
        makeVariantCall("NA00002", new int[] {0}, new double[] {-3, 0}, false, null));
    variantBuilder12.addCalls(
        makeVariantCall("NA00003", new int[] {1}, new double[] {0, -41}, false, null));
    truth.add(variantBuilder12.build());

    return truth;
  }

  public static List<VariantCall> createSomeCalls(List<String> callsetNames, int... genotypes) {
    checkArgument(genotypes.length % 2 == 0);
    checkArgument(callsetNames.size() * 2 >= genotypes.length);
    List<VariantCall> toReturn = new ArrayList<>();
    for (int i = 0; i < genotypes.length; i += 2) {
      toReturn.add(
        VariantCall.newBuilder()
            .setCallSetName(callsetNames.get(i / 2))
            .addGenotype(genotypes[i])
            .addGenotype(genotypes[i + 1])
            .build());
    }

    return toReturn;
  }
}
