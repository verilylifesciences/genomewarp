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
import java.nio.file.Paths;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Helper class for extracting a canonical VCF and a BED file from a gVCF file
 */
public class GvcfToVcfAndBed {

  private static final Logger logger = Logger.getLogger(VcfToVariant.class.getName());

  /**
   * Solve Issue #2 - Extract variant VCF and BED from (tentative) gVCF file
   * 
   * @param rawQueryVcf user-provided (g?)VCF path
   * @param vcfFromGvcf intermediate canonical VCF path
   * @param bedFromGvcf intermediate BED path
   * @return true if query VCF is a gVCF, false otherwise
   */
  public static boolean saveVcfAndBedFromGvcf(String rawQueryVcf, String vcfFromGvcf,
      String bedFromGvcf) throws IOException {
    PrintWriter outBed = new PrintWriter(Files.newBufferedWriter(Paths.get(bedFromGvcf), UTF_8));
    PrintWriter outVcf = new PrintWriter(Files.newBufferedWriter(Paths.get(vcfFromGvcf), UTF_8));
    BufferedReader br = Files.newBufferedReader(Paths.get(rawQueryVcf), UTF_8);
    String line;
    boolean isGvcf = false;
    Long lineNum = 0L;
    while ((line = br.readLine()) != null) {
      // parse line and write it to VCF or BED or none
      lineNum++;
      if (line.length() > 0) {
        if (line.startsWith("#")) {
          // header
          outVcf.println(line);
        } else {
          // not a header - choose destination according to columns
          String[] fields = line.split("\\t");
          if (fields.length < 8) {
            throw new IOException("input VCF file has less than 8 columns");
          }
          String f5 = fields[4];
          if (!".".equals(f5) && !"<*>".equals(f5)) {
            // variant VCF file
            outVcf.println(line);
          }
          if ("PASS".equals(fields[6])) {
            // BED file - compute line
            try {
              Long start = Long.parseLong(fields[1]) - 1L;
              Long end = start + fields[3].length();
              String f8 = fields[7];
              if (!".".equals(f8)) {
                String[] f8Flds = f8.split(";");
                int n = f8Flds.length;
                for (int i = 0; i < n; i++) {
                  String f8Fi = f8Flds[i];
                  if (f8Fi.startsWith("END=")) {
                    String endStr = f8Fi.substring(4);
                    end = Long.parseLong(endStr);
                  }
                }
              }
              StringBuilder sb =
                  new StringBuilder(fields[0]).append("\t").append(start).append("\t").append(end);
              outBed.println(sb.toString());
              // query file is a gVCF
              if (!isGvcf) {
                isGvcf = true;
              }
            } catch (NumberFormatException ex) {
              throw new IOException("malformed integer column in input VCF, line #" + lineNum);
            }
          }
        }
      }
    }
    outBed.close();
    outVcf.close();
    br.close();
    if (!isGvcf) {
      logger.log(Level.INFO, "Query file is not a gVCF - removing intermediate files");
      // remove intermediate files
      Files.delete(Paths.get(vcfFromGvcf));
      // should be empty
      Files.delete(Paths.get(bedFromGvcf));
    }
    return isGvcf;
  }
}
