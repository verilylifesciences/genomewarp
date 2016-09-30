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

import static java.nio.charset.StandardCharsets.UTF_8;

import com.google.common.annotations.VisibleForTesting;
import com.google.genomics.v1.Variant;
import com.google.genomics.v1.VariantCall;
import com.google.protobuf.ListValue;
import com.google.protobuf.Value;

import htsjdk.tribble.TribbleException;
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.IntGenotypeFieldAccessors;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * VcfToVariant contains the method convertVCFtoVariant which takes a
 * VCf file and returns a list of Variant protos. This is done by using
 * htsjdk's {@link VCFFileReader} to parse the VCF then doing field conversion.
 */
public class VcfToVariant {

  private static final Logger logger = Logger.getLogger(VcfToVariant.class.getName());

  private static final IntGenotypeFieldAccessors GENOTYPE_FIELD_ACCESSORS =
      new IntGenotypeFieldAccessors();

  @VisibleForTesting
  static List<String> getAltBases(VariantContext vc) {
    List<Allele> alleleList = vc.getAlternateAlleles();
    List<String> toReturn = new ArrayList<>();

    for (Allele currAllele : alleleList) {
      toReturn.add(currAllele.getBaseString());
    }

    return toReturn;
  }

  @VisibleForTesting
  static List<String> getFilter(VariantContext vc) {
    if (vc.isFiltered()) {
      return (List<String>) ParsingUtils.sortList(vc.getFilters());
    }
    if (vc.filtersWereApplied()) {
      return Arrays.asList(VCFConstants.PASSES_FILTERS_v4);
    }

    return Arrays.asList(VCFConstants.UNFILTERED);
  }

  @VisibleForTesting
  static Map<String, ListValue> getInfo(VariantContext vc, VCFHeader header) {
    Map<String, ListValue> toReturn = new HashMap<>();

    for (Map.Entry<String, Object> entry : vc.getAttributes().entrySet()) {
      String currKey = entry.getKey();
      VCFInfoHeaderLine metaData = header.getInfoHeaderLine(currKey);

      // All info fields must have a corresponding header field.
      if (metaData == null) {
        logger.log(Level.WARNING, String.format("Could not find matching VCF header field, "
            + "skipping info field %s", currKey));
        continue;
      }

      Object currObject = entry.getValue();
      ListValue.Builder listValueBuilder = ListValue.newBuilder();

      VCFHeaderLineType type = metaData.getType();
      if (!(currObject instanceof List)) {
        toReturn.put(currKey,
            listValueBuilder.addValues(createTypedValue(type, currObject)).build());
        continue;
      }

      List<Object> currObjectList = (List<Object>) currObject;
      for (Object currObj : currObjectList) {
        listValueBuilder.addValues(createTypedValue(type, currObj));
      }
      toReturn.put(currKey, listValueBuilder.build());
    }

    return toReturn;
  }

  private static Value createTypedValue(VCFHeaderLineType type, Object value) {
    if (type == VCFHeaderLineType.Flag) {
      return Value.newBuilder().setBoolValue((Boolean) value).build();
    }

    // Booleans are given as Boolean objects. Strangely, Floats and Integers
    // are given as String objects by HTSJDK.
    if (!(value instanceof String)) {
      throw new IllegalStateException("Received non-Boolean, non-List type in non-String format. "
          + "This is most likely due to a change in htsjdk's library.");
    }

    String stringValue = (String) value;
    boolean isNumeric = stringValue.matches("[-+]?\\d+(\\.\\d+)?");

    if (type == VCFHeaderLineType.Integer && isNumeric) {
      return Value.newBuilder().setNumberValue(Integer.parseInt(stringValue)).build();
    }

    if (type == VCFHeaderLineType.Float && isNumeric) {
      return Value.newBuilder().setNumberValue(Double.parseDouble(stringValue)).build();
    }


    return Value.newBuilder().setStringValue(stringValue).build();
  }

  private static Map<Allele, Integer> buildAlleleMap(VariantContext vc) {
    final Map<Allele, Integer> alleleMap = new HashMap<>(vc.getAlleles().size() + 1);
    alleleMap.put(Allele.NO_CALL, -1);

    int alleleNum = 0;
    for (Allele currAllele : vc.getAlleles()) {
      alleleMap.put(currAllele, alleleNum++);
    }

    return alleleMap;
  }

  // Requires non-null accessor
  // Return true if we need to add the listvalue builder to the variant call.
  private static boolean parseInlineGenotypeFields(String field, VariantCall.Builder vcBuilder,
      ListValue.Builder lvBuilder, IntGenotypeFieldAccessors.Accessor accessor, Genotype g) {

    final int[] intValues = accessor.getValues(g);
    if (intValues == null || intValues.length == 0) {
      return false;
    }

    if (field.equals(VCFConstants.GENOTYPE_PL_KEY)) {
      // HTSJDK folds GL's into PL's. We only use PL's to store genotype likelihood.
      for (int i = 0; i < intValues.length; i++) {
        // We add 0.0 to remove the possiblity of getting -0.0.
        vcBuilder.addGenotypeLikelihood(-(double) intValues[i] / 10.0 + 0.0);
      }
      return false;
    }

    for (int i = 0; i < intValues.length; i++) {
      lvBuilder.addValues(Value.newBuilder().setNumberValue(intValues[i]));
    }
    return true;
  }

  // Return true if we need to add the listvalue builder to the variant call.
  private static boolean parseOtherGenotypeFields(String field, VariantContext vc,
      ListValue.Builder lvBuilder, Genotype g, VCFHeader header) {

    if (!g.hasAnyAttribute(field)) {
      return false;
    }

    final VCFFormatHeaderLine metaData = header.getFormatHeaderLine(field);
    if (metaData == null) {
      logger.log(Level.WARNING, String.format("Could not find matching VCF header field for "
          + "genotype field %s", field));
      return false;
    }

    VCFHeaderLineType type = metaData.getType();
    Object value = g.getExtendedAttribute(field);
    final int fieldCount = metaData.getCount(vc);
    if (fieldCount == 1) {
      lvBuilder.addValues(createTypedValue(type, value));
      return true;
    }

    if (!(value instanceof String)) {
      throw new IllegalStateException("received non-Flag genotype field as non-String type");
    }
    String[] valueArray = ((String) value).split(",");
    if (valueArray.length == 1) {
      throw new IllegalStateException(String.format("header indicating a count greater than 1 "
          + "with non-List type found for field %s",
          field));
    }

    boolean allFalse = true;
    for (int i = 0; i < valueArray.length; i++) {
      VCFHeaderLineType thisType = VCFHeaderLineType.String;
      if (!valueArray[i].equals(VCFConstants.MISSING_VALUE_v4)) {
        thisType = type;
        allFalse = false;
      }

      lvBuilder.addValues(createTypedValue(thisType, valueArray[i]));
    }
    // We only add the lvBuilder if there is at least one non-missing value
    return !allFalse;
  }

  @VisibleForTesting
  static List<VariantCall> getCalls(VariantContext vc, VCFHeader header) {
    List<VariantCall> toReturn = new ArrayList<>();

    for (String currSample : header.getGenotypeSamples()) {
      if (!vc.hasGenotype(currSample)) {
        continue;
      }

      Genotype currGenotype = vc.getGenotype(currSample);
      VariantCall.Builder vcBuilder = VariantCall.newBuilder();
      vcBuilder.setCallSetName(currSample);

      // Get GT info.
      final Map<Allele, Integer> alleleStrings = buildAlleleMap(vc);
      vcBuilder.addGenotype(alleleStrings.get(currGenotype.getAllele(0)));
      for (int i = 1; i < currGenotype.getPloidy(); i++) {
        vcBuilder.addGenotype(alleleStrings.get(currGenotype.getAllele(i)));
      }

      // Set phasing (not applicable to haploid).
      if (currGenotype.isPhased() && currGenotype.getPloidy() > 1) {
        vcBuilder.setPhaseset("*");
      }

      // Get rest of the genotype info.
      Map<String, ListValue> genotypeInfo = new HashMap<>();

      // Set filters
      if (currGenotype.isFiltered()) {
        genotypeInfo.put(VCFConstants.GENOTYPE_FILTER_KEY, ListValue.newBuilder()
            .addValues(Value.newBuilder().setStringValue(currGenotype.getFilters()).build())
            .build());
      }

      for (final String field : vc.calcVCFGenotypeKeys(header)) {
        // We've already handled genotype
        if (field.equals(VCFConstants.GENOTYPE_KEY)) {
          continue;
        }

        ListValue.Builder listValueBuilder = ListValue.newBuilder();
        if (field.equals(VCFConstants.GENOTYPE_FILTER_KEY)) {
          // This field has already been dealt with
          continue;
        } else {
          final IntGenotypeFieldAccessors.Accessor accessor =
              GENOTYPE_FIELD_ACCESSORS.getAccessor(field);

          if (accessor != null) {
            // The field is a default inline field.
            if (!parseInlineGenotypeFields(field, vcBuilder, listValueBuilder, accessor,
                currGenotype)) {
              continue;
            }
          } else {
            // Other field, we'll get type/other info from header.
            if (!parseOtherGenotypeFields(field, vc, listValueBuilder, currGenotype,
                header)) {
              continue;
            }
          }
        }

        genotypeInfo.put(field, listValueBuilder.build());
      }

      vcBuilder.putAllInfo(genotypeInfo);
      toReturn.add(vcBuilder.build());
    }

    return toReturn;
  }

  private static boolean validVersion(File filepath) throws IOException {
    BufferedReader reader = Files.newBufferedReader(filepath.toPath(), UTF_8);

    // The first line must be the header
    String firstLine = reader.readLine();
    reader.close();

    try {
      VCFHeaderVersion version = VCFHeaderVersion.getHeaderVersion(firstLine);

      // If the version is greater than or equal to 4.2, we cannot handle it
      if (version.isAtLeastAsRecentAs(VCFHeaderVersion.VCF4_2)) {
        return false;
      }
    } catch (TribbleException.InvalidHeader msg) {
      throw new IOException(msg);
    }

    return true;
  }

  public static List<Variant> convertContextToVariant(Iterable<VariantContext> variantContexts,
      VCFHeader header) {
    List<Variant> finalList = new ArrayList<>();

    for (final VariantContext vc : variantContexts) {
      String referenceName = vc.getChr();

      // htsjdk start/end are both 1-based closed
      long variantStart = vc.getStart() - 1;
      long variantEnd = vc.getEnd();

      String refBases = vc.getReference().getBaseString();
      boolean hasQuality = vc.hasLog10PError();
      double quality = vc.getPhredScaledQual();

      boolean hasId = vc.hasID();
      String[] ids = vc.getID().split(",");

      List<String> altBases = getAltBases(vc);

      List<String> filters = getFilter(vc);

      Map<String, ListValue> info = getInfo(vc, header);

      List<VariantCall> calls = getCalls(vc, header);

      // Actually build the variant
      Variant.Builder variantBuilder = Variant.newBuilder().setReferenceName(referenceName).
          setStart(variantStart).setEnd(variantEnd).setReferenceBases(refBases);
      if (hasQuality) {
        variantBuilder.setQuality(quality);
      }
      if (hasId) {
        for (String id : ids) {
          variantBuilder.addNames(id);
        }
      }
      for (String base : altBases) {
        variantBuilder.addAlternateBases(base);
      }
      for (VariantCall call : calls) {
        variantBuilder.addCalls(call);
      }
      for (String filter : filters) {
        variantBuilder.addFilter(filter);
      }
      variantBuilder.putAllInfo(info);

      finalList.add(variantBuilder.build());
    }

    return finalList;
  }


  /**
   * Converts a VCF at the given file path to a list of Variant protos.
   * <p> This function reads from a path to a VCF file using the
   * HTSJDK library. Following that, it performs field conversion
   * to create a Variant from the resulting {@link VCFFileReader} object.
   * Note that this implementation depends heavily on HTSJDK and makes the
   * same assumptions as HTSJDK (e.g. integer values GQ). Furthermore, due
   * to using HTSJDK, this tool is limited to working with VCF versions of
   * at most 4.1, and will return null if given a VCF version higher than 4.1.
   *
   * @param filepath Path to the VCF file
   * @return A list of Variant protos
   */
  public static List<Variant> convertVcfToVariant(File filepath) throws IOException {
    return convertVcfToVariant(filepath, true);
  }

  public static List<Variant> convertVcfToVariant(File filepath, boolean checkVersion)
      throws IOException {
    // Check version
    if (checkVersion && !validVersion(filepath)) {
      return null;
    }

    VCFFileReader vcfReader = new VCFFileReader(filepath, false);
    VCFHeader header = vcfReader.getFileHeader();

    return convertContextToVariant(vcfReader, header);
  }

}
