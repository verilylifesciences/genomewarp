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
import static org.junit.Assert.fail;

import com.google.genomics.v1.Variant;
import com.google.genomics.v1.VariantCall;
import com.google.protobuf.ListValue;
import com.verily.genomewarp.utils.VcfToVariant;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import java.io.File;
import java.util.List;
import java.util.Map;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
//@MediumTest(MediumTestAttribute.FILE)
public final class VcfToVariantTest {

  private static final String VALID_VCF_4_1 = "utils/valid-4.1.vcf";

  private static final String EXTERNAL_VCF = "utils/smallNA12877_S1.genome.vcf";

  private static final boolean SUPPORTS_V4_2 = false;

  private static final List<Variant> TRUTH = getTruth(SUPPORTS_V4_2);

  @Test
  public void testGetAltBases() throws Exception {
    File vcfFile =
        new File(VcfToVariant.class.getClassLoader().getResource(VALID_VCF_4_1).getFile());
    VCFFileReader vcfReader = new VCFFileReader(vcfFile, false);
    int currVariant = 0;

    for (final VariantContext vc : vcfReader) {
      List<String> altBases = VcfToVariant.getAltBases(vc);
      assertEquals(altBases, TRUTH.get(currVariant).getAlternateBasesList());
      currVariant++;
    }
  }

  @Test
  public void testGetFilter() throws Exception {
    File vcfFile =
        new File(VcfToVariant.class.getClassLoader().getResource(VALID_VCF_4_1).getFile());
    VCFFileReader vcfReader = new VCFFileReader(vcfFile, false);
    int currVariant = 0;

    for (final VariantContext vc : vcfReader) {
      List<String> filters = VcfToVariant.getFilter(vc);
      assertEquals(filters, TRUTH.get(currVariant).getFilterList());
      currVariant++;
    }
  }

  @Test
  public void testGetInfo() throws Exception {
    File vcfFile =
        new File(VcfToVariant.class.getClassLoader().getResource(VALID_VCF_4_1).getFile());
    VCFFileReader vcfReader = new VCFFileReader(vcfFile, false);
    VCFHeader vcfHeader = vcfReader.getFileHeader();
    int currVariant = 0;

    for (final VariantContext vc : vcfReader) {
      Map<String, ListValue> info = VcfToVariant.getInfo(vc, vcfHeader);
      assertEquals(info, TRUTH.get(currVariant).getInfo());
      currVariant++;
    }
  }

  @Test
  /**
   * Note: HTSJDK cannot distinguish between PASS filters and
   * unfiltered for variant calls (it can for general filters).
   * As a result, all variant calls which PASS do not have filter
   * information applied.
   */
  public void testGetCalls() throws Exception {
    File vcfFile =
        new File(VcfToVariant.class.getClassLoader().getResource(VALID_VCF_4_1).getFile());
    VCFFileReader vcfReader = new VCFFileReader(vcfFile, false);
    VCFHeader vcfHeader = vcfReader.getFileHeader();
    int currVariant = 0;

    for (final VariantContext vc : vcfReader) {
      List<VariantCall> callList = VcfToVariant.getCalls(vc, vcfHeader);
      assertEquals(callList, TRUTH.get(currVariant).getCallsList());
      currVariant++;
    }
  }

  @Test
  public void testCompleteRun() throws Exception {
    File vcfFile =
        new File(VcfToVariant.class.getClassLoader().getResource(VALID_VCF_4_1).getFile());
    assertEquals(TRUTH, VcfToVariant.convertVcfToVariant(vcfFile));
  }

  @Test
  public void testSuccessfulExternalVCF() {
    File vcfFile =
        new File(VcfToVariant.class.getClassLoader().getResource(EXTERNAL_VCF).getFile());
    try {
      VcfToVariant.convertVcfToVariant(vcfFile);
    } catch (Exception ex) {
      fail("conversion threw exception when expected to run successfully");
    }
  }
}
