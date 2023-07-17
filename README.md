# tbHCA
HCA (**H**ybrid **C**apture **A**nalysis) is a pipeline written in Nextflow that analyzes hybridization capture sequencing data. 

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

# Benchling 

Make sure you have your Benchling API credentials in a file named .env

```
WAREHOUSE_USERNAME = 'username'
WAREHOUSE_PASSWORD = 'password'
WAREHOUSE_URL = 'postgres-warehouse.tome.benchling.com'

API_KEY = 'api_key_here'
API_URL = 'https://tome.benchling.com/'
```

If it's your first time attempting to query Benchling via Python, you might run into an issue regarding the `root.crt` file. Please reference this [docs link](https://help.benchling.com/hc/en-us/articles/9714802961421-Access-your-data-warehouse) to fix this issue. 
Never publicly upload your API credentials! 

# Dependencies

The following section lists the dependencies for this project.

## System Requirements

- Python 3.10

## Packages

Make sure you have the following packages installed:

- [minimap2](https://github.com/lh3/minimap2)
- [bwa](https://github.com/lh3/bwa)
- [CRISPResso2](https://github.com/pinellolab/CRISPResso2)
- [samtools](https://github.com/samtools/samtools)
- [pysam](https://pypi.org/project/pysam/)
- [biopython](https://pypi.org/project/biopython/)
- [glob2](https://pypi.org/project/glob2/)
- [pandas](https://pypi.org/project/pandas/)
- [psycopg2](https://pypi.org/project/psycopg2/)
- [nextflow](https://github.com/nextflow-io/nextflow)
- [fastp](https://github.com/OpenGene/fastp)
- [dotenv](https://github.com/motdotla/dotenv)

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





