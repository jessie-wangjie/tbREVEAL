# tbHCA
HCA (**H**ybrid **C**apture **A**nalysis) is a pipeline written in Nextflow that analyzes hybridization capture sequencing data. 

# Dependencies

The following section lists the dependencies for this project.

## System Requirements

- Python 3.10

## Packages

Make sure you have the following packages installed:

- [minimap2 2.26](https://github.com/lh3/minimap2)
- [bwa 0.7.17](https://github.com/lh3/bwa)
- [CRISPResso2 2.2.12](https://github.com/pinellolab/CRISPResso2)
- [samtools 1.17](https://github.com/samtools/samtools)
- [pysam 0.21.0](https://pypi.org/project/pysam/)
- [biopython 1.81](https://pypi.org/project/biopython/)
- [glob2 0.7](https://pypi.org/project/glob2/)
- [pandas 2.0.3](https://pypi.org/project/pandas/)
- [psycopg2 2.9.6](https://pypi.org/project/psycopg2/)
- [nextflow 23.04.1](https://github.com/nextflow-io/nextflow)
- [fastp 0.23.4](https://github.com/OpenGene/fastp)
- [dotenv 1.0.0](https://github.com/motdotla/dotenv)

Note that the versions listed are the versions when I made the pipeline, but they will probably run in future versions as well (unless there's a major change to one of them, causing dependency issues.) 

All of these can be installed via conda (or mamba, as I prefer):

```
conda create -n HCA -c bioconda -c conda-forge -c anaconda -c bioconda minimap2 bwa crispresso2 samtools pysam biopython glob2 pandas psycopg2 nextflow fastp python-dotenv
```

OR 


```
mamba create -n HCA -c bioconda -c conda-forge -c anaconda -c bioconda minimap2 bwa crispresso2 samtools pysam biopython glob2 pandas psycopg2 nextflow fastp python-dotenv
```

then

```
conda activate HCA
```

# Quick start

After cloning the repo, you can run the pipeline using:

```
nextflow run pipeline.nf \
--r1 R1.fastq.gz \
--r2 R2.fastq.gz \
--metadata metadata.csv \
--reference hg38.fa \
--outdir my_output_directory \
--sample_name my_sample_name
```

Make sure to use absolute paths when referencing a filename. A full example may look like this:

```
nextflow run pipeline.nf \
--r1 /data/reads/JM_Custom_TB000174c_20230613/S4-up-350-rep1_L1_ds.1258ae405c464bd1b100534fdd5e2699/S4-up-350-rep1_S1_L001_R1_001.fastq.gz \
--r2 /data/reads/JM_Custom_TB000174c_20230613/S4-up-350-rep1_L1_ds.1258ae405c464bd1b100534fdd5e2699/S4-up-350-rep1_S1_L001_R2_001.fastq.gz \
--metadata /data/CM_Custom_TB000175c_20230613_analysis/20230329-final_panel_up.csv \
--reference /data/references/hg38.fa \
--outdir /data/JM_Custom_TB000174c_20230613_analysis/ \
--sample_name S4-up-350-rep1
```

## Metadata file

The metadata file should be in .csv format and looks like:

| Chromosome | Start    | Stop     | Name     | Score | Strand | Code | Seq                                                           | len | GC   | mtDNA/rDNA hit | Target | Target chr | Target Start | Target Stop | Gap Size |
|------------|----------|----------|----------|-------|--------|------|---------------------------------------------------------------|-----|------|----------------|--------|------------|--------------|-------------|----------|
| 15         | 44711415 | 44711535 | AA1520-UP| 1     | +      | good | CACTGCGTCGCTGGCTTGGAGACAGGTGACGGTCCCTGCGGGCCTTGTCCTGATTGGCTGGGCACGCGTTTAATATAAGTGGAGGCGTCGCGCTGGCGGGCATTCCTGAAGCTGACAGCA | 120 | 61.67 | no             | AA1520 | NA         | NA           | NA          | NA       |
| 16         | 10895583 | 10895703 | AA1522-UP| 1     | +      | good | CCAGCCCTGCCCCGCCTCTCCCTCGTTCCCCACCAGCCCTCTTTCCAGAAATTTCCTTCTTCATCCAAGGGACTTTTCCTCCCAGAACCCGACACAGACACCATCAACTGCGACCAGTTC | 120 | 57.5  | no             | AA1522 | NA         | NA           | NA          | NA       |
| 1          | 26317940 | 26318060 | AA1542-UP| 1     | +      | good | AACCAAAAGAAGCCTCCAGACAGCCCTGAGATCACCTAAAAAGCTGCTACCAAGACAGCCACGAAGATCCTACCAAAATGAAGCGCTTCCTCTTCCTCCTACTCACCATCAGCCTCCTGG | 120 | 50.83 | no             | AA1542 | NA         | NA           | NA          | NA       |
| 3          | 50611519 | 50611639 | AA1544-UP| 1     | +      | good | AGTCACCTCTGGCCCGTCAAGCCCTCCCAATGCCCGGCAGCTAGCACGAAGCCCCTGTTCTCCCGTGCGCCCCTCGTGGTGGCCGGGAAGGGGGCAGAGAGCCGCGCTTACCCCTGAACG | 120 | 69.17 | no             | AA1544 | NA         | NA           | NA          | NA       |
| X          | 139536457| 139536577| AA729-UP | 1     | +      | good | CCTATAACACTTGCCAACCAAAGGTGCTGTTGATCTGAAATTGCTTTTTTAAATTAATGCAGTGATTTTTCTTTAACATCTAGTGACAGACACTGGGGTCACATTTGCAGCTGGACCATA | 120 | 38.33 | no             | AA729  | NA         | NA           | NA          | NA       |
| 1          | 184150972| 184151092| CAS002-UP| 1     | -      | good | CTTACACTACTTGCTTCAATGACTTTGAACTTGGCGTGCCGTTCTGTGGCTTTGCTGCTGTCTGAATCACATGCTTTTGCCTGCATTACCAAGCAGGGCTTGGAGCCCAGTCTCAGGAGG | 120 | 50.83 | no             | CAS002 | chr1       | 184151093    | 184151093   | 1        |

The most important columns are: Target, Chromosome, Start, Stop, Strand. The pipeline should work with just these columns.  

## Benchling 

Make sure you have your Benchling API credentials in a file named .env inside the bin/utils/ folder

```
WAREHOUSE_USERNAME = 'username'
WAREHOUSE_PASSWORD = 'password'
WAREHOUSE_URL = 'postgres-warehouse.tome.benchling.com'

API_KEY = 'api_key_here'
API_URL = 'https://tome.benchling.com/'
```

If it's your first time attempting to query Benchling via Python, you might run into an issue regarding the `root.crt` file. Please reference this [docs link](https://help.benchling.com/hc/en-us/articles/9714802961421-Access-your-data-warehouse) to fix this issue. 
Never publicly upload your API credentials! 

## Output

The pipeline creates symlinks to the work directory to the following files/directories:

`alignments`
SAM files containing each site-specific read alignments

`amplicons`
Sequences used for each site (attL, attR, cryptic B)

`attL_extracted_reads`
For sites where attL was detected, output a truncated version of the read to match attL. Used for input into CRISPResso2 

`attR_extracted_reads` 
For sites where attR was detected, output a truncated version of the read to match attL. Used for input into CRISPResso2 

`cs2_output`
CRISPResso2 output

`extracted_reads`
FASTQ files for every site. This is whats used to create the "alignments" directory

`indel_info`
Information on indels

`initial_alignment` 
Initial alignment of all reads to the reference genome

`integration_and_indel_stats.csv`
The most important output -- a csv file containing each site and a bunch of important information including integration percentage, indel percentage, risk category, gene, number of reads... etc

`integration_stats.csv`
A .csv file similar to above, except it just contains information on integration percentage. This will likely be removed in a future version. 

`qc`
Output from fastp containing QC info

`target_info`
A .csv containing information on each site including the cryptic sequences



