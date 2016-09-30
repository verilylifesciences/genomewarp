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

import static org.junit.Assert.assertEquals;

import com.verily.genomewarp.utils.Fasta;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.experimental.runners.Enclosed;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(Enclosed.class)
public class FastaTest {

  private static final String FASTA_PATH = "utils/test.fa";

  private static final String DOS_FASTA_PATH = "utils/testDOS.fa";

  @RunWith(JUnit4.class)
  public static final class LoadFilesTest {

    private static String chr2Bases, chr3Bases, chr17Bases, chr20Bases;

    @BeforeClass
    public static void setup() {
      LoadFilesTest.chr2Bases = "AAGGATGACAAGAAATAAATAGCTACCGCTATTATGAGTCGCATGGTAAG";
      LoadFilesTest.chr3Bases = "AACGTTTNNX";
      LoadFilesTest.chr17Bases = "GCAGGGGCCACGGGGGGAGCAGCCTCTGGCATTCTGGGAGCTTCATCTGG"
          + "ACCTGGGTCTTCAGTGAACCATTGTTCAATATCGTCCGGGGACAGCATCA"
          + "AATCATCCAT";
      LoadFilesTest.chr20Bases = "TGtgC";
    }

    @Test
    public void testLoadSequence() throws IOException {
      Map<String, String> expected = new HashMap<String, String>();
      expected.put("chr2", LoadFilesTest.chr2Bases.toUpperCase());
      expected.put("chr3", LoadFilesTest.chr3Bases.toUpperCase());
      expected.put("chr17", LoadFilesTest.chr17Bases.toUpperCase());
      expected.put("chr20", LoadFilesTest.chr20Bases.toUpperCase());

      Fasta fastaTest =
          new Fasta(FastaTest.class.getClassLoader().getResource(FASTA_PATH).getFile());
      assertEquals(fastaTest.get("chr2", -1, -1), LoadFilesTest.chr2Bases.toUpperCase());
      assertEquals(fastaTest.get("chr3", -1, -1), LoadFilesTest.chr3Bases.toUpperCase());
      assertEquals(fastaTest.get("chr17", -1, -1), LoadFilesTest.chr17Bases.toUpperCase());
      assertEquals(fastaTest.get("chr20", -1, -1), LoadFilesTest.chr20Bases.toUpperCase());
    }
  }

  @RunWith(JUnit4.class)
  public static final class RangeGetTest {

    @Test
    public void testGet() throws IOException {
      Fasta fastaTest =
          new Fasta(FastaTest.class.getClassLoader().getResource(FASTA_PATH).getFile());
      assertEquals(fastaTest.get("chr2", 3, 10), "GATGACA");
      assertEquals(fastaTest.get("chr2", 3, 4), "G");
      assertEquals(fastaTest.get("chr2", -1, 10), "AAGGATGACA");
      assertEquals(fastaTest.get("chr2", 35, -1), "GAGTCGCATGGTAAG");
      assertEquals(fastaTest.get("chr17", 0, 6), "GCAGGG");
      assertEquals(fastaTest.get("chr17", 100, 199), "AATCATCCAT");
      assertEquals(fastaTest.get("chr20", -1, -1), "TGTGC");
      assertEquals(fastaTest.get("should be empty", 10, 23), "");
    }

    @Test(expected = IllegalArgumentException.class)
    public void testDosEnding() throws IOException {
      Fasta fastaTest =
          new Fasta(FastaTest.class.getClassLoader().getResource(DOS_FASTA_PATH).getFile());
      fastaTest.get("chr17", 3, 10);
    }
  }
}
