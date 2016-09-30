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

import com.google.genomics.v1.Variant;
import com.google.genomics.v1.VariantCall;
import com.google.protobuf.ListValue;
import com.google.protobuf.Value;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;

import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * VariantToVcf contain the method convertVariantToVcf which takes a Variant
 * proto and writes a VCF file to the given path. This is done by first
 * constructing VariantContexts then using htsjdk's VariantContextWriter
 * to write the header and each line of the VCF file
 */
public class VariantToVcf {

  private static final String SOURCE = "variant_to_vcf";

  private static VariantContextWriter vcfWriter = null;

  private static String getKind(Value.KindCase kind) {
    switch (kind) {
        case NULL_VALUE: return "NULL_VALUE";
        case NUMBER_VALUE: return "NUMBER_VALUE";
        case STRING_VALUE: return "STRING_VALUE";
        case BOOL_VALUE: return "BOOL_VALUE";
        case STRUCT_VALUE: return "STRUCT_VALUE";
        case LIST_VALUE: return "LIST_VALUE";
        default: return "KIND_NOT_SET";
    }
  }

  private static Object getNumberValue(double numValue, VCFHeaderLineType type) {
    if (String.valueOf(numValue).matches("[-+]?\\d+(\\.\\d+)?")) {
      switch (type) {
          case Integer: return String.format("%d", (int) numValue);
          case Float: return String.format("%s", numValue);
          default: break;
      }
    }
    return String.valueOf(numValue);
  }

  private static Object parseObject(String key, Value in, VCFHeaderLineType type) {
    // Case on type
    switch (in.getKindCase()) {
        case NULL_VALUE: throw new IllegalStateException(String.format("field %s contained "
            + "a null value", key));
        case NUMBER_VALUE: return getNumberValue(in.getNumberValue(), type);
        case STRING_VALUE: return in.getStringValue();
        case BOOL_VALUE: return (Boolean) in.getBoolValue();
        default: throw new IllegalStateException(String.format("field %s contained a %s type, which"
            + " is not supported by VariantToVcf", key, getKind(in.getKindCase())));
    }
  }

  private static int[] glsToPls(List<Double> gls) {
    int[] pls = new int[gls.size()];
    Double adjust = Double.NEGATIVE_INFINITY;
    for (Double gl : gls) {
      adjust = Math.max(adjust, gl);
    }

    int i = 0;
    for (Double gl : gls) {
      pls[i++] = (int) Math.round(Math.min(-10 * (gl - adjust), Integer.MAX_VALUE));
    }

    return pls;
  }

  private static Genotype getGenotypeFromVariantCall(VCFHeader header, VariantCall in,
      List<Allele> refList) {
    List<Allele> list = new ArrayList<>();
    for (int allele : in.getGenotypeList()) {
      if (allele == -1) {
        list.add(Allele.NO_CALL);
        continue;
      }

      list.add(refList.get(allele));
    }
    GenotypeBuilder builder = new GenotypeBuilder(in.getCallSetName(), list);

    // Set phasing
    String phaseSet = in.getPhaseset();
    if (phaseSet != null && !phaseSet.isEmpty()) {
      builder.phased(true);
    }

    // Set PL using
    int numGL = in.getGenotypeLikelihoodCount();
    if (numGL == 0) {
      builder.noPL();
    } else {
      builder.PL(glsToPls(in.getGenotypeLikelihoodList()));
    }

    // Set other attributes
    Map<String, Object> attributes = new HashMap<>();
    for (Map.Entry<String, ListValue> entry : in.getInfo().entrySet()) {
      String currKey = entry.getKey();
      List<Value> currValue = entry.getValue().getValuesList();
      int currSize = currValue.size();
      if (currSize == 0) {
        // If the key is present in the info map, there must
        // be non-zero attributes associated with it
        throw new IllegalStateException(String.format("genotype field %s has length 0", currKey));
      }

      // Handle the in line fields, which are all integers
      if (currKey.equals("GQ") || currKey.equals("AD") || currKey.equals("DP")) {
        if (currValue.get(0).getKindCase() != Value.KindCase.NUMBER_VALUE) {
          throw new IllegalStateException(String.format("genotype field %s must have a number value"
              + ", had %s", currKey, getKind(currValue.get(0).getKindCase())));
        }

        // GQ and DP should be length 1
        if (currKey.equals("GQ") || currKey.equals("DP")) {
          if (currSize != 1) {
            throw new IllegalStateException(String.format("genotype field %s should have length 1, "
                + "had length %d", currKey, currValue.size()));
          }
        }

        // Convert all values into integers
        int[] intList = new int[currSize];
        int i = 0;
        for (Value val : currValue) {
          intList[i++] = ((Double) val.getNumberValue()).intValue();
        }

        if (currKey.equals("GQ")) {
          builder.GQ(intList[0]);
        } else if (currKey.equals("DP")) {
          builder.DP(intList[0]);
        } else if (currKey.equals("AD")) {
          builder.AD(intList);
        }

        continue;
      } else if (currKey.equals(VCFConstants.GENOTYPE_FILTER_KEY)) {
        List<String> filters = new ArrayList<>();
        for (Value val : currValue) {
          filters.add(val.getStringValue());
        }
        builder.filters(filters);
        continue;
      }

      // Handle extraneous fields
      VCFHeaderLineType type = header.getFormatHeaderLine(currKey).getType();
      if (currSize == 1) {
        attributes.put(currKey, parseObject(currKey, currValue.get(0), type));
      } else {
        List<Object> objectList = new ArrayList<>();
        for (Value val : currValue) {
          // Handles where only some list values are missing
          if (val.getKindCase() == Value.KindCase.STRING_VALUE) {
            if (val.getStringValue().equals(VCFConstants.MISSING_VALUE_v4)) {
              objectList.add(VCFConstants.MISSING_VALUE_v4);
              continue;
            }
          }
          objectList.add(parseObject(currKey, val, type));
        }
        attributes.put(currKey, objectList);
      }
    }

    if (attributes.size() == 0) {
      builder.noAttributes();
    } else {
      builder.attributes(attributes);
    }

    return builder.make();
  }

  private static VariantContext getContextFromVariant(VCFHeader header, Variant in) {
    // Create allele list
    String refBases = in.getReferenceBases();
    List<String> altBases = in.getAlternateBasesList();
    List<Allele> alleleList = new ArrayList<>();
    alleleList.add(Allele.create(refBases.getBytes(StandardCharsets.UTF_8), true));
    for (String base : altBases) {
      alleleList.add(Allele.create(base.getBytes(StandardCharsets.UTF_8), false));
    }

    // Note htsjdk start/end are both 1-based closed
    VariantContextBuilder builder = new VariantContextBuilder(SOURCE, in.getReferenceName(),
        in.getStart() + 1, in.getEnd(), alleleList);

    // Set Quality
    if (in.getQuality() != 0.0) {
      builder.log10PError(-in.getQuality() / 10);
    }

    // Set names
    int numNames = in.getNamesCount();
    int i = 0;
    String namesString = "";
    for (String name : in.getNamesList()) {
      if (i == numNames - 1) {
        break;
      }
      namesString += name + ",";
    }
    if (numNames == 0) {
      builder.noID();
    } else {
      // Also handle fence post
      builder.id(namesString + in.getNames(numNames - 1));
    }

    // Set filters
    boolean hadFilter = false;
    for (String filter : in.getFilterList()) {
      hadFilter = true;
      if (filter.equals(VCFConstants.PASSES_FILTERS_v4)) {
        builder.passFilters();
        break;
      }
      if (filter.equals(VCFConstants.MISSING_VALUE_v4)) {
        builder.unfiltered();
        break;
      }

      builder.filter(filter);
    }
    if (!hadFilter) {
      builder.unfiltered();
    }

    // Set info
    Map<String, Object> toSet = new HashMap<>();
    for (Map.Entry<String, ListValue> entry : in.getInfo().entrySet()) {
      String currKey = entry.getKey();
      List<Value> currValue = entry.getValue().getValuesList();
      int currSize = currValue.size();

      if (currSize == 0) {
        // If the key is present in the info map, there must
        // be non-zero attributes associated with it
        throw new IllegalStateException(String.format("info field %s has length 0", currKey));
      }

      VCFHeaderLineType type = header.getInfoHeaderLine(currKey).getType();
      if (currSize == 1) {
        toSet.put(currKey, parseObject(currKey, currValue.get(0), type));
        continue;
      }

      List<Object> objectList = new ArrayList<>();
      for (Value val : currValue) {
        objectList.add(parseObject(currKey, val, type));
      }
      toSet.put(currKey, objectList);
    }
    builder.attributes(toSet);

    // Set calls
    List<Genotype> genotypes = new ArrayList<>();
    for (VariantCall vc : in.getCallsList()) {
      genotypes.add(getGenotypeFromVariantCall(header, vc, alleleList));
    }
    if (genotypes.size() == 0) {
      builder.noGenotypes();
    } else {
      builder.genotypes(genotypes);
    }

    return builder.make();
  }

  /**
   * Converts a list of variants into a VCF file, given a {@link VCFHeader}
   * and an {@link OutputStream}.
   * <p> This function uses HTSJDK to create {@link VariantCall} objects then
   * writes them to the given output stream. Note that this implementation depends
   * heavily on HTSJDK and makes the same assumptions as HTSJDK (e.g. integer values GQ).
   *
   * @param header The header to use to generate the output VCF
   * @param variants A list of variants to encode in the output VCF
   * @param os The output stream to which to write the generated VCF
   */
  public static void convertVariantToVcf(VCFHeader header, List<Variant> variants,
      OutputStream os, boolean writeHeader) {
    if (vcfWriter == null) {
      vcfWriter = new VariantContextWriterBuilder().clearOptions()
          .setOutputVCFStream(os).build();
    }

    if (writeHeader) {
      vcfWriter.writeHeader(header);
    }

    for (Variant currVariant : variants) {
      vcfWriter.add(getContextFromVariant(header, currVariant));
    }
  }
}
