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

import static com.verily.genomewarp.utils.GenomeWarpTestUtils.makeVariant;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.genomics.v1.Variant;
import com.verily.genomewarp.ReferenceGenomeDifference;
import java.util.Arrays;
import org.junit.Rule;
import org.junit.Test;
import org.junit.experimental.runners.Enclosed;
import org.junit.rules.ExpectedException;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameter;
import org.junit.runners.Parameterized.Parameters;

@RunWith(Enclosed.class)
public final class ReferenceGenomeDifferenceTest {

  @RunWith(Parameterized.class)
  public static final class ValidCreateTests {
    @Parameters(name = "{index}: create({0})={1}")
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            {1L, "ACGT", "A"}, {3L, "T", "C"}, {4L, "A", "AT"},
          });
    }

    @Parameter(0)
    public long expectedPosition;

    @Parameter(1)
    public String expectedQuery;

    @Parameter(2)
    public String expectedTarget;

    @Test
    public void createReferenceGenomeDifference() {
      ReferenceGenomeDifference actualRgd =
          ReferenceGenomeDifference.create(expectedPosition, expectedQuery, expectedTarget);
      assertTrue(actualRgd.queryPosition() == expectedPosition);
      assertTrue(actualRgd.queryDna().equals(expectedQuery));
      assertTrue(actualRgd.targetDna().equals(expectedTarget));
    }
  }

  @RunWith(Parameterized.class)
  public static final class InvalidCreateTests {
    @Parameters(name = "{index}: create({0})")
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            {1L, "A", "A"},
            {3L, "", "C"},
            {4L, "C", ""},
            {5L, "", ""},
            {6L, "A", "TCA"},
            {7L, "AC", "TG"},
            {8L, "ACC", "ATG"},
            {9L, "AT", "ATC"},
          });
    }

    @Parameter(0)
    public long position;

    @Parameter(1)
    public String query;

    @Parameter(2)
    public String target;

    @Rule public final ExpectedException thrown = ExpectedException.none();

    @Test
    public void createInvalidReferenceGenomeDifference() {
      thrown.expect(IllegalArgumentException.class);
      ReferenceGenomeDifference.create(position, query, target);
    }
  }

  @RunWith(Parameterized.class)
  public static final class ClassificationTests {
    @Parameters(name = "{index}: classify({0})")
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            {1L, "A", "C", true, false, false},
            {2L, "A", "ACCC", false, true, false},
            {3L, "AC", "A", false, false, true},
            {4L, "TGGC", "T", false, false, true},
          });
    }

    @Parameter(0)
    public long position;

    @Parameter(1)
    public String query;

    @Parameter(2)
    public String target;

    @Parameter(3)
    public boolean expectedIsSnv;

    @Parameter(4)
    public boolean expectedIsInsertion;

    @Parameter(5)
    public boolean expectedIsDeletion;

    @Test
    public void classifyReferenceGenomeDifference() {
      ReferenceGenomeDifference actualRgd =
          ReferenceGenomeDifference.create(position, query, target);
      assertTrue(actualRgd.isSnv() == expectedIsSnv);
      assertTrue(actualRgd.isInsertion() == expectedIsInsertion);
      assertTrue(actualRgd.isDeletion() == expectedIsDeletion);
    }
  }

  public static final class ClassifyNoDifferenceTests {
    @Test
    public void classifyNoDifference() {
      assertFalse(ReferenceGenomeDifference.NO_DIFFERENCE.isSnv());
      assertFalse(ReferenceGenomeDifference.NO_DIFFERENCE.isInsertion());
      assertFalse(ReferenceGenomeDifference.NO_DIFFERENCE.isDeletion());
    }
  }

  @RunWith(Parameterized.class)
  public static final class OverlapTests {
    @Parameters(name = "{index}: overlapsVariant({0}, {1})={2}")
    public static Iterable<Object[]> testData() {
      return Arrays.asList(
          new Object[][] {
            {2L, "A", "C", makeVariant("chr1", 1, "A", "C"), false},
            {2L, "A", "C", makeVariant("chr1", 2, "A", "C"), true},
            {2L, "A", "C", makeVariant("chr1", 3, "A", "C"), false},
            {2L, "A", "C", makeVariant("chr1", 1, "A", "ATT"), false},
            {2L, "A", "C", makeVariant("chr1", 2, "A", "ATT"), true},
            {2L, "A", "C", makeVariant("chr1", 3, "A", "ATT"), false},
            {5L, "A", "C", makeVariant("chr1", 1, "ACAG", "A"), false},
            {5L, "A", "C", makeVariant("chr1", 2, "ACAG", "A"), true},
            {5L, "A", "C", makeVariant("chr1", 3, "ACAG", "A"), true},
            {5L, "A", "C", makeVariant("chr1", 4, "ACAG", "A"), true},
            {5L, "A", "C", makeVariant("chr1", 5, "ACAG", "A"), true},
            {5L, "A", "C", makeVariant("chr1", 6, "ACAG", "A"), false},
            {8L, "A", "ACACA", makeVariant("chr1", 7, "A", "C"), false},
            {8L, "A", "ACACA", makeVariant("chr1", 8, "A", "C"), true},
            {8L, "A", "ACACA", makeVariant("chr1", 9, "A", "C"), false},
            {11L, "ATG", "A", makeVariant("chr1", 10, "A", "C"), false},
            {11L, "ATG", "A", makeVariant("chr1", 11, "A", "C"), true},
            {11L, "ATG", "A", makeVariant("chr1", 12, "A", "C"), true},
            {11L, "ATG", "A", makeVariant("chr2", 13, "A", "C"), true},
            {11L, "ATG", "A", makeVariant("chr1", 14, "A", "C"), false},
          });
    }

    @Parameter(0)
    public long position;

    @Parameter(1)
    public String query;

    @Parameter(2)
    public String target;

    @Parameter(3)
    public Variant variant;

    @Parameter(4)
    public boolean expectedOverlaps;

    @Test
    public void referenceGenomeDifferenceOverlapsVariant() {
      ReferenceGenomeDifference diff = ReferenceGenomeDifference.create(position, query, target);
      assertTrue(diff.overlapsVariant(variant) == expectedOverlaps);
    }
  }
}
