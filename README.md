# tbHCA
REVEAL is a pipeline written in Nextflow that analyzes sequencing data (such as Hybrid Capture and LMPCR) which captures multiple types of edits in one run.

# Dependencies

The easiest way to get all dependencies needed is by using Docker.

First, clone the repo.

Second, install Docker on your instance.

Once docker is installed, you can build the Docker image by first navigating to the top of the code repository, then:

```docker build -t hca:latest .```

If you wish to directly use a package manager like conda or mamba, there is an environment.yml file that you can utilize.

# Quick start

Your nextflow.config file should look something like this. Most important things are the project id, benchling credentials, aws credentials, and basespace credentials.

```
report.enabled = true
report.overwrite = true
docker.enabled = true
process.container = 'hca:latest'

params {
    manifest {
        name = "HCA"
        description = "This pipeline analyzes hybridization capture sequencing data for both on-and-off target recombination events, and outputs key information"
        author = "Thomas Biondi (thomas.biondi@tome.bio)"
        version = "1.0.0"
    }

    project_id='<project id on benchling, i.e. CTB029_1>'
    outdir = "output"
    initial_mapper = 'bwa'

    bucket_name='s3://tb-ngs-genomics-quilt/'

    BENCHLING_WAREHOUSE_USERNAME='<benchling warehouse username>'
    BENCHLING_WAREHOUSE_PASSWORD='<benchling warehouse password>'
    BENCHLING_WAREHOUSE_URL='postgres-warehouse.tome.benchling.com'
    BENCHLING_API_URL='https://tome.benchling.com/'
    BENCHLING_API_KEY='<benchling api key>'

    BS_ACCESS_TOKEN = '<basespace access token>'
    BS_API_SERVER = 'https://api.basespace.illumina.com'

    AWS_ACCESS_KEY_ID = '<aws access key id>'
    AWS_SECRET_ACCESS_KEY = '<aws secret access key>'

}
```


