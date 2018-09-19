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

import com.google.common.collect.ComparisonChain;

/**
 * class GenomeRange stores genomic ranges for a specified chromosome. Furthermore
 * it also stores information for liftover, such as the score and strand. GenomeRange
 * implements Comparable, to allow sorting of GenomeRanges.
 */
public class GenomeRange implements Comparable<GenomeRange>{
  private String chromosome;
  private long start;
  private long end;
  private String name;
  // negative score indicates no score
  private int score;
  // true -> positive strand, false -> negative strand
  private boolean isPositiveStrand;

  public GenomeRange(String chromosome, long start, long end) {
    this.chromosome = chromosome;
    this.start = start;
    this.end = end;
    this.name = "";
    this.score = -1;
    this.isPositiveStrand = true;
  }

  public GenomeRange(String chromosome, long start, long end, String name, boolean strand) {
    this.chromosome = chromosome;
    this.start = start;
    this.end = end;
    this.name = name;
    this.score = -1;
    this.isPositiveStrand = strand;
  }

  public String getChromosome() {
    return chromosome;
  }

  public long getStart() {
    return start;
  }

  public long getEnd() {
    return end;
  }

  public String getName() {
    return name;
  }

  public void setName(String name) {
    this.name = name;
  }

  public int getScore() {
    return score;
  }

  public boolean isPositiveStrand() {
    return isPositiveStrand;
  }

  @Override
  public int compareTo(GenomeRange range) {
    return ComparisonChain.start()
        .compare(chromosome, range.chromosome)
        .compare(start, range.start)
        .compare(end, range.end)
        .compare(isPositiveStrand, range.isPositiveStrand())
        .result();
  }

  public boolean overlaps(GenomeRange range) {
    return range != null && chromosome.equals(range.getChromosome()) && end > range.getStart()
        && start < range.getEnd();
  }

  public boolean isSameInterval(GenomeRange range) {
    return range != null && chromosome.equals(range.getChromosome()) && (start == range.getStart())
        && (end == range.getEnd());
  }

  public boolean includes(GenomeRange range) {
    return range != null && chromosome.equals(range.getChromosome()) && start <= range.getStart()
        && end >= range.getEnd();
  }

  public GenomeRange getIntersection(GenomeRange range) {
    if (!this.overlaps(range)) {
      return null;
    }
    return new GenomeRange(chromosome, Math.max(start, range.getStart()),
        Math.min(end, range.getEnd()));
  }

  @Override
  public String toString() {
    return String.format("%s\t%s\t%s\t%s\t.\t%s", chromosome, start, end, name,
        isPositiveStrand ? "+" : "-");
  }
}
