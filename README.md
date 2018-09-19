## Disclaimer

This is not an official Verily product.

# GenomeWarp

GenomeWarp is a command-line tool that translates genetic variants in
confidently-called genomic regions from one genome assembly version to another,
such as from GRCh37 to GRCh38.

## Purpose

The goal of GenomeWarp is to translate the variation within a set of regions
deemed "confidently-called" in one genome assembly to another genome assembly.
In cases where a VCF file represents "all variation in an individual with
respect to the genome assembly against which the VCF was generated", GenomeWarp
can be used to transform that data to the analogous set of all variation in the
individual with respect to a new genome assembly.

This is semantically different from existing tools that support the translation
of VCF files from one assembly to another, including:

*   [NCBI Remap](http://www.ncbi.nlm.nih.gov/genome/tools/remap)
*   [CrossMap](http://crossmap.sourceforge.net/)

These tools operate only on the sites present in an input VCF, and return the
representation of those sites in a new genome assembly. This does not capture
all variation, however. Consider an individual who has sequence reads that
indicate they match the GRCh37 reference genome assembly at position
[GRCh37.chr1:169,519,049](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr1%3A169519049-169519049)
(i.e. the individual's genotype is T/T). Because the individual is homozygous
reference at that site, there will be no variation present in their VCF file
created on GRCh37. However, the analogous position on the updated GRCh38
reference genome assembly, position
[GRCh38.chr1:169,549,811](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr1%3A169549811-169549811),
has the reference base C. Consequently, if the individual's read data were
analyzed on GRCh38, they would be identified as homozygous for a C->T SNP.
Because this site is not present in the input GRCh37 VCF, it is never added
when creating a GRCh38 VCF by these other tools.

## Nomenclature and format definitions

In the below descriptions, the genome assembly on which the confidently-called
regions and variants are given is denoted the "query" assembly. The genome
assembly onto which the user wishes to warp the variants and regions is denoted
the "target" assembly.

File formats:

*   [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
*   [Chain](https://genome.ucsc.edu/goldenPath/help/chain.html)
*   [FASTA](https://en.wikipedia.org/wiki/FASTA_format)
*   [VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf)

## Inputs

To warp variants from a query assembly to a target assembly, five inputs are
required:

*   A BED file of confidently-called regions in the query assembly
*   A VCF file containing containing variants in the query assembly
    *   Note: This may include variants outside of the confidently-called
        regions, but those variants will be ignored.
*   A FASTA file containing the query assembly sequence
*   A FASTA file containing the target assembly sequence
*   A Chain file providing coordinate transformation from query to target
    assembly

The BED file of confidently-called regions can be created by emitting an
all-sites VCF file (when calling variants with the
[GATK](https://software.broadinstitute.org/gatk/)) and filtering the
homozygous-reference calls at a desired quality threshold.

## Common file downloads
FASTA files can be downloaded from NCBI or the UCSC Genome Browser, among other
places. See the [NCBI How
To](http://www.ncbi.nlm.nih.gov/guide/howto/dwn-genome/) for details.

Common chain file downloads are available from the UCSC Genome Browser at the
URLs

```
http://hgdownload.cse.ucsc.edu/goldenpath/${QUERY}/liftOver/${QUERY}To${TARGET^}.over.chain.gz
```

for the appropriate definitions of `QUERY` and `TARGET`. For example, the chain
to transform from hg19 to hg38 is
[http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz).

## Building this project

1. git clone this repository.
1. Use a recent version of [Apache Maven](http://maven.apache.org/download.cgi)
(e.g., version 3.3.3) to build this code:

```
mvn package
```

## Running the GenomeWarp tool

Once all five input files are available, performing the transformation involves
running a single Java program. The driver script is `GenomeWarpSerial.java`, and
it generates two output files:

* A BED file of all confidently-called regions in the target genome assembly
* A VCF file containing all variants in the confidently-called regions in the
  target assembly.

The program is executed as follows:

```bash
java -jar target/verilylifesciences-genomewarp-1.1.0-runnable.jar \
  --lift_over_chain_path "${chain}" \
  --raw_query_vcf "${queryvcf}" \
  --raw_query_bed "${querybed}" \
  --ref_query_fasta "${queryfasta}" \
  --ref_target_fasta "${targetfasta}" \
  --work_dir "${workdir}" \
  --output_variants_file "${targetvcf}" \
  --output_regions_file "${targetbed}"
```

When run, logging statements provide progress indications. GenomeWarp should
convert a single-sample VCF containing millions of variants genome-wide in under
30 minutes.

## Notes
There are multiple reasons why a confidently-called region in the query assembly
(and any variants therein) may not appear in the target assembly. GenomeWarp is
deliberately conservative in tricky cases, preferring to omit a
confidently-called region and its constituent variants if there is not an
unambiguous mapping. The guarantee GenomeWarp provides is that all
confidently-called regions in the target assembly faithfully reproduce the same
haplotypes as were provided in the query assembly (i.e., GenomeWarp gives 100%
specificity at a possible sacrifice to sensitivity).

GenomeWarp currently handles variant-only VCF files (i.e. gVCFs are not
supported). A gVCF can be processed using the workaround [described
here](https://github.com/verilylifesciences/genomewarp/issues/2).
