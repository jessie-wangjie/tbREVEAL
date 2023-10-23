# tbHCA
HCA (**H**ybrid **C**apture **A**nalysis) is a pipeline written in Nextflow that analyzes hybridization capture sequencing data. 

# Dependencies

The easiest way to get all dependencies needed is by using Docker. 

First, clone the repo. 

Second, install Docker on your instance. 

Once docker is installed, you can build the Docker image by first navigating to the top of the code repository, then: 

```docker build -t hca:latest .```

If you wish to directly use a package manager like conda or mamba, there is an environment.yml file that you can utilize. 

# Quick start

There are multiple ways to run the pipeline. The easiest way is to directly edit the nextflow.config file. By default, Nextflow will use this file to find the input files. 

```
nextflow run pipeline.nf -with-docker hca
```

Alternatively, you can give Nextflow the input files at runtime as such (this will overwrite any parameters given in the nextflow.config file):

```
nextflow run pipeline.nf \
--samplesheet samplesheet.csv \
--reference hg38.fa \
--outdir my_output_directory \
--deduplication_method LOOSE \
--collapse_condition CARGO \
--project_name my_project_name
```

Right now, use relative paths when specifying inputs. I am currently working to add support for both absolute and relative paths. 

## Samplesheet

HCA utilizes a samplesheet, which contains information about all samples within a particular project. An example samplesheet is in the examples folder. 

The samplesheet should have the following columns: sample_name,fastq_dir,probe_list,group

sample_name is just the name of your sample, can be anything you want

fastq_dir is the relative path to your FASTQ files. 

probe_list is the relative path to the probe information file. More information about that can be seen below. 

group is specifying whether the sample is a control, treated sample, or anything else. 

## Probe information file

One of the key files is the probe information sheet. This file should be in .csv format and looks like:

| Chromosome | Start    | Stop     | Name     | Score | Strand | Code | Seq                                                           | len | GC   | mtDNA/rDNA hit | Target | Target chr | Target Start | Target Stop | Gap Size |
|------------|----------|----------|----------|-------|--------|------|---------------------------------------------------------------|-----|------|----------------|--------|------------|--------------|-------------|----------|
| 15         | 44711415 | 44711535 | AA1520-UP| 1     | +      | good | CACTGCGTCGCTGGCTTGGAGACAGGTGACGGTCCCTGCGGGCCTTGTCCTGATTGGCTGGGCACGCGTTTAATATAAGTGGAGGCGTCGCGCTGGCGGGCATTCCTGAAGCTGACAGCA | 120 | 61.67 | no             | AA1520 | NA         | NA           | NA          | NA       |
| 16         | 10895583 | 10895703 | AA1522-UP| 1     | +      | good | CCAGCCCTGCCCCGCCTCTCCCTCGTTCCCCACCAGCCCTCTTTCCAGAAATTTCCTTCTTCATCCAAGGGACTTTTCCTCCCAGAACCCGACACAGACACCATCAACTGCGACCAGTTC | 120 | 57.5  | no             | AA1522 | NA         | NA           | NA          | NA       |
| 1          | 26317940 | 26318060 | AA1542-UP| 1     | +      | good | AACCAAAAGAAGCCTCCAGACAGCCCTGAGATCACCTAAAAAGCTGCTACCAAGACAGCCACGAAGATCCTACCAAAATGAAGCGCTTCCTCTTCCTCCTACTCACCATCAGCCTCCTGG | 120 | 50.83 | no             | AA1542 | NA         | NA           | NA          | NA       |
| 3          | 50611519 | 50611639 | AA1544-UP| 1     | +      | good | AGTCACCTCTGGCCCGTCAAGCCCTCCCAATGCCCGGCAGCTAGCACGAAGCCCCTGTTCTCCCGTGCGCCCCTCGTGGTGGCCGGGAAGGGGGCAGAGAGCCGCGCTTACCCCTGAACG | 120 | 69.17 | no             | AA1544 | NA         | NA           | NA          | NA       |
| X          | 139536457| 139536577| AA729-UP | 1     | +      | good | CCTATAACACTTGCCAACCAAAGGTGCTGTTGATCTGAAATTGCTTTTTTAAATTAATGCAGTGATTTTTCTTTAACATCTAGTGACAGACACTGGGGTCACATTTGCAGCTGGACCATA | 120 | 38.33 | no             | AA729  | NA         | NA           | NA          | NA       |
| 1          | 184150972| 184151092| CAS002-UP| 1     | -      | good | CTTACACTACTTGCTTCAATGACTTTGAACTTGGCGTGCCGTTCTGTGGCTTTGCTGCTGTCTGAATCACATGCTTTTGCCTGCATTACCAAGCAGGGCTTGGAGCCCAGTCTCAGGAGG | 120 | 50.83 | no             | CAS002 | chr1       | 184151093    | 184151093   | 1        |

The most important columns are: Target, Chromosome, Start, Stop, Strand. The pipeline should work with just these columns.  

An example probe information file is in the examples folder.  

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

under construction



