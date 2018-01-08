/*
 * Copyright 2016 Verily Life Sciences LLC.
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
 * in compliance with the License. You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software distributed under the License
 * is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
 * or implied. See the License for the specific language governing permissions and limitations under
 * the License.
 */

package com.verily.genomewarp.utils;

import static java.nio.charset.StandardCharsets.UTF_8;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Helper class for extracting a canonical VCF and a BED file from a gVCF file.
 */
public class GvcfToVcfAndBed {
  // class constants
  private static final Logger logger = Logger.getLogger(GvcfToVcfAndBed.class.getName());
  
  private static final int MIN_VCF_FIELDS = 8;
  private static final int CHROM_FIELD_INDEX = 0;
  private static final int START_FIELD_INDEX = 1;
  private static final int END_FIELD_INDEX = 3;
  private static final int ALT_BASES_FIELD_INDEX = 4;
  private static final int FILTER_FIELD_INDEX = 6;
  private static final int FORMAT_FIELD_INDEX = 7;
  private static final String UNKNOWN_VALUE = ".";
  private static final String GVCF_ALT_ALLELE = "<*>";
  private static final String HEADER_START_STR = "#";
  private static final String FIELD_SEP_REGEX = "\\t";
  private static final String PASS_FILTER_VALUE = "PASS";
  private static final String FORMAT_FIELD_SEP_REGEX = ";";
  private static final String FORMAT_END_KEY = "END";
  private static final String BED_FIELD_SEP_STR = "\t";
  
  /**
   * Solve Issue #2 - Extract variant VCF and BED from (tentative) gVCF file
   * 
   * @param rawQueryGvcf user-provided gVCF path
   * @param vcfFromGvcf intermediate canonical VCF path
   * @param bedFromGvcf intermediate BED path
   * @return true if query gVCF was parsed successfully, false otherwise
   */
  public static boolean saveVcfAndBedFromGvcf(String rawQueryGvcf, String vcfFromGvcf,
      String bedFromGvcf) {
    String line;
    boolean successfullyParsedGvcf = true;
    Long lineNum = 0L;
    Path bedFromGvcfPath = null;
    Path vcfFromGvcfPath = null;
    Path rawQueryGvcfPath = null;
    try {
      bedFromGvcfPath = Paths.get(bedFromGvcf);
      vcfFromGvcfPath = Paths.get(vcfFromGvcf);
      rawQueryGvcfPath = Paths.get(rawQueryGvcf);
    } catch (Exception ex) {
      // at least one invalid path
      logger.log(Level.SEVERE, "Caught exception: " + ex.toString());
      return false;
    }

    try (PrintWriter outVcf = new PrintWriter(Files.newBufferedWriter(vcfFromGvcfPath, UTF_8));
        PrintWriter outBed = new PrintWriter(Files.newBufferedWriter(bedFromGvcfPath, UTF_8));
        BufferedReader br = Files.newBufferedReader(rawQueryGvcfPath, UTF_8)) {
      while ((line = br.readLine()) != null) {
        // parse line and write it to VCF or BED or none
        lineNum++;
        if (line.startsWith(HEADER_START_STR)) {
          // header
          outVcf.println(line);
        } else {
          // not a header - choose destination according to columns
          String[] fields = line.split(FIELD_SEP_REGEX);
          if (fields.length < MIN_VCF_FIELDS) {
            throw new IOException("input VCF file has less than " + MIN_VCF_FIELDS + " columns,"
                + " line #" + lineNum);
          }
          String altBases = fields[ALT_BASES_FIELD_INDEX];
          if (isVariantRecord(altBases)) {
            // variant VCF file
            outVcf.println(line);
          }
          if (PASS_FILTER_VALUE.equals(fields[FILTER_FIELD_INDEX])) {
            // BED file - compute line
            try {
              outBed.println(bedLineFromVcfFields(fields));
            } catch (NumberFormatException ex) {
              throw new IOException("malformed integer column in input VCF, line #" + lineNum);
            }
          }
        }
      }
    } catch (Exception ex) {
      logger.log(Level.SEVERE, "Caught exception: " + ex.toString());
      successfullyParsedGvcf = false;
    }
    if (!successfullyParsedGvcf) {
      logger.log(Level.INFO,
          "gVCF query file processing error or not a gVCF - removing intermediate files");
      // remove intermediate files
      if (Files.exists(vcfFromGvcfPath)) {
        try {
          Files.delete(vcfFromGvcfPath);
        } catch (IOException ex) {
          logger.log(Level.INFO, "Failed to delete VCF intermediate file: " + ex.toString());
        }
      }
      if (Files.exists(bedFromGvcfPath)) {
        try {
          Files.delete(bedFromGvcfPath);
        } catch (IOException ex) {
          logger.log(Level.INFO, "Failed to delete BED intermediate file: " + ex.toString());
        }
      }
    }
    return successfullyParsedGvcf;
  }

  /**
   * @param altBases alternate bases field value
   * @return true if the field value denotes a strictly variant record
   */
  private static boolean isVariantRecord(String altBases) {
    return !UNKNOWN_VALUE.equals(altBases) && !GVCF_ALT_ALLELE.equals(altBases);
  }

  /**
   * @param fields the gVCF fields
   * @return the corresponding BED line
   */
  private static String bedLineFromVcfFields(String[] fields) throws NumberFormatException {
    // VCF is 1-based full closed indexing, whereas BED is 0-based half open indexing,
    // hence the -1L here.
    Long start = Long.parseLong(fields[START_FIELD_INDEX]) - 1L;
    Long end = start + fields[END_FIELD_INDEX].length();
    String rawFormatField = fields[FORMAT_FIELD_INDEX];
    String[] formatFieldTokens = rawFormatField.split(FORMAT_FIELD_SEP_REGEX);
    int n = formatFieldTokens.length;
    for (int i = 0; i < n; i++) {
      String token = formatFieldTokens[i];
      if (token.startsWith(FORMAT_END_KEY + "=")) {
        String endStr = token.substring(FORMAT_END_KEY.length() + 1);
        end = Long.parseLong(endStr);
      }
    }
    StringBuilder sb = new StringBuilder(fields[CHROM_FIELD_INDEX]).append(BED_FIELD_SEP_STR)
        .append(start).append(BED_FIELD_SEP_STR).append(end);
    return sb.toString();
  }
}
