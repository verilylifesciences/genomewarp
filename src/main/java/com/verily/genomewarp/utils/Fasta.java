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

import com.google.common.collect.ImmutableMap;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Class for interacting with sequences from a single genome assembly FASTA.
 */
public class Fasta implements AutoCloseable {

  private static Logger logger = Logger.getLogger(Fasta.class.getName());

  private String fastaName;

  private Map<String, String> data;

  private final Map<String, Long> jumpTable;

  private Map<String, Long> chrSize;

  public static final String MISSING_CHROMOSOME = "";

  public Fasta(String fastaName) throws IOException {
    this.fastaName = fastaName;
    this.data = new HashMap<>();
    this.jumpTable = createJumpTable();
  }

  // Creates a jump table from chromosome to location in the file
  private Map<String, Long> createJumpTable() throws IOException {
    ImmutableMap.Builder<String, Long> builder = ImmutableMap.builder();
    chrSize = new HashMap<>();
    BufferedReader thisReader =
        new BufferedReader(new InputStreamReader(new FileInputStream(this.fastaName), UTF_8));

    String line;
    Long position = 0L;
    Long length = 0L;
    String lastChr = "";
    while ((line = thisReader.readLine()) != null) {
      // We add 1 to account for new line characters
      if (line.length() > 0 && line.charAt(0) == '>') {
        String chromosome = line.split("\\s+")[0].replaceAll("^>", "");
        builder.put(chromosome, position);

        if (!lastChr.isEmpty()) {
          chrSize.put(lastChr, length);
          length = 0L;
        }
        lastChr = chromosome;
      } else {
        length += line.length();
      }
      position += line.length() + 1;
    }

    // Fencepost
    if (!lastChr.isEmpty()) {
      chrSize.put(lastChr, length);
    }

    return builder.build();
  }

  // For batch loading
  private Map<String, String> loadSequenceFromFile(List<String> toLoad) throws IOException {
    String line;

    Map<String, String> dataMap = new HashMap<>();

    for (String chromosome : toLoad) {
      // Get position from jump table
      if (!this.jumpTable.containsKey(chromosome)) {
        continue;
      }
      Long jumpTo = this.jumpTable.get(chromosome);

      try (BufferedReader thisReader =
          new BufferedReader(new InputStreamReader(new FileInputStream(this.fastaName), UTF_8))) {
        StringBuilder pieces = new StringBuilder();
        thisReader.skip(jumpTo);
        // The first line must be the chromosome header
        line = thisReader.readLine();
        if (line == null || !chromosome.equals(line.split("\\s+")[0].replaceAll("^>", ""))) {
          throw new IllegalArgumentException("Unsupported fasta file: newline encoding must be '\r'"
              + " or '\n'");
        }

        while ((line = thisReader.readLine()) != null) {
          if (line.length() > 0 && line.charAt(0) == '>') {
            dataMap.put(chromosome, pieces.toString());
            break;
          } else {
            pieces.append(line.trim().toUpperCase());
          }
        }
        // Fencepost
        if (!chromosome.equals("")) {
          dataMap.put(chromosome, pieces.toString());
        }
      }
    }

    return dataMap;
  }

  public void preload(List<String> toLoad) throws IOException {
    Map<String, String> newData = new HashMap<>();
    List<String> finalLoad;
    if (data.size() > 0) {
      // If we have some data already, only load what we have to
      finalLoad = new ArrayList<>();
      Set<String> exists = data.keySet();
      for (String currChr : toLoad) {
        if (exists.contains(currChr)) {
          newData.put(currChr, this.data.get(currChr));
        } else {
          // Not present in our map already, need to load
          finalLoad.add(currChr);
        }
      }
      if (finalLoad.size() == 0) {
        return;
      }
    } else {
      finalLoad = toLoad;
    }

    newData.putAll(loadSequenceFromFile(finalLoad));
    this.data = newData;
  }

  /**
   * Returns the string bases associated with a specified chromsome in a given range.
   * <p> This function returns a string of bases which correspond to the given
   * chromosome between the start (inclusive) and end (exclusive) indices. Furthermore,
   * either bound can be set to -1 to ignore that bound. If the range or chromosome
   * is invalid, this function returns the empty string.
   *
   * @param chromosome The desired chromosome
   * @param start The inclusive, 0-based start index
   * @param end The exclusive, 0-based end index
   * @return a string containing the bases in the [start, end) range
   */
  public String get(String chromosome, long start, long end) {
    String chrString = MISSING_CHROMOSOME;
    try {
      if (this.jumpTable.containsKey(chromosome) && !this.data.containsKey(chromosome)) {
        preload(Arrays.asList(chromosome));
      }
      chrString = this.data.get(chromosome);
    } catch (IOException ex) {
      throw new RuntimeException("failed to read from specified FASTA file: "
          + ex.getMessage());
    }

    if (chrString == null) {
      logger.log(Level.WARNING, "Failed to load chromosome " + chromosome);
      // Store the missing chromosome with "" to prevent rereading entire file
      data.put(chromosome, MISSING_CHROMOSOME);
      return MISSING_CHROMOSOME;
    }

    // If we get a cached missing value, immediately return it
    if (MISSING_CHROMOSOME.equals(chrString)) {
      return MISSING_CHROMOSOME;
    }

    if (start == -1) {
      start = 0;
    }
    if (end == -1 || end > chrString.length()) {
      end = chrString.length();
    }

    if (start > Integer.MAX_VALUE || end > Integer.MAX_VALUE) {
      throw new IllegalArgumentException("Fasta.get(): canot specify bases indices greater than "
          + "integer max values");
    }

    return chrString.substring((int) start, (int) end);
  }

  /**
   * Returns the length of the contig associated with the given chromosome
   * or -1 if the specified chromosome does not exist in the fasta.
   *
   * @param chromosome The desired chromosome
   * @return The length of the contig
   */
  public long getChromosomeSize(String chromosome) {
    if (!chrSize.containsKey(chromosome)) {
      return -1;
    }

    return chrSize.get(chromosome);
  }

  /**
   * Returns the reference names in the FASTA file loaded into the
   * jump table by their ordering in the original FASTA.
   *
   * @return A list of reference names by original FASTA position
   */
  public Map<String, Long> getReferenceNameAndLengthInInputOrder() {
    List<Map.Entry<String, Long>> entryList = new ArrayList<>(jumpTable.entrySet());
    Collections.sort(entryList, new Comparator<Object>() {
      @SuppressWarnings("unchecked")
      public int compare(Object v1, Object v2) {
        return ((Map.Entry<String, Long>) v1).getValue()
            .compareTo(((Map.Entry<String, Long>) v2).getValue());
      }
    });

    // Use linked hash map to preserve order
    Map<String, Long> toReturn = new LinkedHashMap<>();
    for (int i = 0; i < entryList.size(); i++) {
      Map.Entry<String, Long> currEntry = entryList.get(i);
      toReturn.put(currEntry.getKey(), chrSize.get(currEntry.getKey()));
    }

    return toReturn;
  }

  /**
   * Closes associated resources
   */
  @Override
  public void close() {
      return;
  }
}
