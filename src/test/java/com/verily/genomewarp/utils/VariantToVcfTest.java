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

import static com.verily.genomewarp.utils.GenomeWarpTestUtils.getTruth;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.genomics.v1.Variant;
import com.verily.genomewarp.utils.VariantToVcf;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
public final class VariantToVcfTest {

  private static final String VALID_VCF_HEADER = "utils/valid-4.2-header.vcf";

  private static final boolean SUPPORTS_V4_2 = false;

  private static final List<Variant> TRUTH = getTruth(SUPPORTS_V4_2);

  private static final int BUFFER_SIZE = 1024;

  // Note these are slightly different from the VCF lines in testdata/valid-4.1.vcf in the
  // following ways:
  // 1. Since protobuf has no sense of "ordering" for info/format fields, these fields are
  //    alphabetized in the desired output.
  // 2. No trailing missing values
  // 3. GL's are converted to equivalent PL's
  private static final String[] VCF_LINES = {
    // Note on this line, GQ higher than 99 (in this case 200) is rounded to 99
    "19	14370	rs6054257	G	A	29	PASS	AF=0.5;DB;DP=14;H2;NS=3	GT:DP:GQ:HQ	0|0:1:99:51,"
        + "51	1|0:8:48:51,51	1/1:5:43",
    "20	17330	.	T	A	3	q10	AF=0.017;DP=11;NS=3	GT:DP:FT:GQ:HQ	0|0:3:q10:49:58,50	0|1:5:s50:3:.,"
        + "3	0/0:3:q10:41",
    "20	1110696	rs6040355	A	G,T	67	PASS	AA=T;AF=0.333,0.667;DB;DP=10;NS=2	GT:DP:GQ:HQ	1|2:6:21:23"
        + ",27	2|1:0:2:18,2	2/2:4:35",
    "20	1230237	.	T	.	47	PASS	AA=T;DP=13;NS=3	GT:DP:GQ:HQ	0|0:7:54:56,60	0|0:4:48:51,51	0/"
        + "0:2:61",
    "20	1234567	microsat1	GTC	G,GTCTC	50	PASS	AA=G;DP=9;NS=3	GT:DP:GQ	0/1:4:35	0/2:2:17	1/1:"
        + "3:40",
    "20	2234568	.	C	TC	50	PASS	AA=G;DP=9;NS=3;SVTYPE=BND	GT:DP:GQ	0/1:4:35	0/1:2:17	1/1:3:40",
    "20	2234569	.	C	CT	50	PASS	AA=G;DP=9;NS=3;SVTYPE=BND	GT:DP:GQ	0/1:4:35	0/1:2:17	./.:3:40",
    "Y	17330	.	T	A	3	q10	AF=0.017;DP=11;NS=3	GT:PL	0:490,0	0:30,0	1:0,410"
  };

  private static final List<String> STRING_HEADER = new ArrayList<String>();

  @Test
  public void testVariantToVcf() throws Exception {
    // Load HEADER for conversions
    File vcfFile = new File(VariantToVcfTest.class.getClassLoader().getResource(VALID_VCF_HEADER)
        .getFile());
    Scanner scan = new Scanner(vcfFile);

    while (scan.hasNextLine()) {
      STRING_HEADER.add(scan.nextLine());
    }

    VCFHeader vcfHeader = new VCFFileReader(vcfFile, false).getFileHeader();
    ByteArrayOutputStream baos = new ByteArrayOutputStream(BUFFER_SIZE);
    int vcfLineNum = 0;

    // Convert Variant
    for (Variant vcfLine : TRUTH) {
      VariantToVcf.convertVariantToVcf(vcfHeader, Arrays.asList(vcfLine), baos, vcfLineNum == 0);

      String[] res = baos.toString(StandardCharsets.UTF_8.toString()).split("\\r?\\n");
      List<String> headerLines = new ArrayList<String>();

      int i = 0;
      int length = res.length;
      for (String line : res) {
        // The first length - 1 lines should be header
        if (i++ == length - 1) {
          break;
        }

        headerLines.add(line);
      }

      // The input header must be a subset of the output header
      assertTrue(headerLines.containsAll(STRING_HEADER));

      // Now compare the actual VCF line
      assertEquals(res[length - 1], VCF_LINES[vcfLineNum++]);
    }
  }
}
