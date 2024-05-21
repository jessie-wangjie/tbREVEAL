nextflow.enable.dsl=2


params.outdir = "output"
params.collapse_condition = 'Complete'
params.initial_mapper = 'bwa'
params.notebook_template = "${workflow.projectDir}/bin/report_generation.ipynb"
params.bam2html_path = "${workflow.projectDir}/bin/utils/bam2html.py"
params.dinucleotides = ''
params.cosmic_info = "/data/cryptic_prediction/data/cosmic/cancer_gene_census.csv"
params.BENCHLING_WAREHOUSE_USERNAME=''
params.BENCHLING_WAREHOUSE_PASSWORD=''
params.BENCHLING_WAREHOUSE_URL=''
params.BENCHLING_WAREHOUSE_API_KEY=''
params.BENCHLING_WAREHOUSE_SDK_KEY=''
params.BENCHLING_API_URL=''
params.bucket_name='s3://tb-ngs-genomics-quilt/'
params.project_id=''
params.quilt_package_name="HybridCapture/${params.project_id}"

process GET_PROJECT_INFO {
    input:
        val project_id
        path raw_reads
        val benchling_warehouse_username
        val benchling_warehouse_password
        val benchling_warehouse_url
        val benchling_sdk_api_key
        val benchling_api_url
    output:
        path('samplesheet.csv')
    script:
        """
        export WAREHOUSE_USERNAME='${benchling_warehouse_username}'
        export WAREHOUSE_PASSWORD='${benchling_warehouse_password}'
        export WAREHOUSE_URL='${benchling_warehouse_url}'
        export API_KEY='${benchling_sdk_api_key}'
        export API_URL='${benchling_api_url}'
        get_project_info.py --project_id ${project_id}
        """
}

process DOWNLOAD_GTEX_DATA {
    input:

    output:
        path('gtex_gene_median_tpm.csv')
    script:
        """
        wget https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
        gunzip GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
        mv GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct gtex_gene_median_tpm.csv
        """
}

process DOWNLOAD_READS {
    input:
        val project_id
    output:
        path "*"
    script:
        """
        run_id=\$(bs list projects --filter-term=^${project_id}\$ -f csv | grep "${project_id}" | cut -f 2 -d ',' | tail -1)
        bs download projects -i \${run_id} -o . --extension=fastq.gz --no-metadata
        mv */* .
        find . -type d -empty -exec rmdir {} +
        """
}

process DOWNLOAD_REFERENCE_GENOME {
    input:
        val reference_species
    output:
        path("*.{fa,fna,fasta}"), emit: reference_fasta
        path("*.{fa,fna,fasta}.*"), emit: reference_index
    script:

    if (reference_species == "Human" || reference_species == "Homo sapiens") {
        species_reference_fasta_path = 's3://tomebfx-data/references/hg38_no_alts/hg38_no_alts.fa'
        species_reference_amb_path = 's3://tomebfx-data/references/hg38_no_alts/hg38_no_alts.fa.amb'
        species_reference_ann_path = 's3://tomebfx-data/references/hg38_no_alts/hg38_no_alts.fa.ann'
        species_reference_bwt_path = 's3://tomebfx-data/references/hg38_no_alts/hg38_no_alts.fa.bwt'
        species_reference_fai_path = 's3://tomebfx-data/references/hg38_no_alts/hg38_no_alts.fa.fai'
        species_reference_pac_path = 's3://tomebfx-data/references/hg38_no_alts/hg38_no_alts.fa.pac'
        species_reference_sa_path = 's3://tomebfx-data/references/hg38_no_alts/hg38_no_alts.fa.sa'
    } else if (reference_species == "Mouse") {
        species_reference_fasta_path = 's3://tomebfx-data/references/hg38_no_alts/hg38_no_alts.fa'
        species_reference_amb_path = 's3://tomebfx-data/references/hg38_no_alts/hg38_no_alts.fa.amb'
        species_reference_ann_path = 's3://tomebfx-data/references/hg38_no_alts/hg38_no_alts.fa.ann'
        species_reference_bwt_path = 's3://tomebfx-data/references/hg38_no_alts/hg38_no_alts.fa.bwt'
        species_reference_fai_path = 's3://tomebfx-data/references/hg38_no_alts/hg38_no_alts.fa.fai'
        species_reference_pac_path = 's3://tomebfx-data/references/hg38_no_alts/hg38_no_alts.fa.pac'
        species_reference_sa_path = 's3://tomebfx-data/references/hg38_no_alts/hg38_no_alts.fa.sa'
    } else if (reference_species == "Monkey" || reference_species == "Macaca fascicularis" || reference_species == "NHP") {
        species_reference_fasta_path = 's3://tomebfx-data/references/Macaca_fascicularis_6/GCA_011100615.1_Macaca_fascicularis_6.0_genomic.fna'
        species_reference_amb_path = 's3://tomebfx-data/references/Macaca_fascicularis_6/GCA_011100615.1_Macaca_fascicularis_6.0_genomic.fna.amb'
        species_reference_ann_path = 's3://tomebfx-data/references/Macaca_fascicularis_6/GCA_011100615.1_Macaca_fascicularis_6.0_genomic.fna.ann'
        species_reference_bwt_path = 's3://tomebfx-data/references/Macaca_fascicularis_6/GCA_011100615.1_Macaca_fascicularis_6.0_genomic.fna.bwt'
        species_reference_fai_path = 's3://tomebfx-data/references/Macaca_fascicularis_6/GCA_011100615.1_Macaca_fascicularis_6.0_genomic.fna.fai'
        species_reference_pac_path = 's3://tomebfx-data/references/Macaca_fascicularis_6/GCA_011100615.1_Macaca_fascicularis_6.0_genomic.fna.pac'
        species_reference_sa_path = 's3://tomebfx-data/references/Macaca_fascicularis_6/GCA_011100615.1_Macaca_fascicularis_6.0_genomic.fna.sa'
    }
    """
    aws s3 cp ${species_reference_fasta_path} .
    aws s3 cp ${species_reference_amb_path} .
    aws s3 cp ${species_reference_ann_path} .
    aws s3 cp ${species_reference_bwt_path} .
    aws s3 cp ${species_reference_fai_path} .
    aws s3 cp ${species_reference_pac_path} .
    aws s3 cp ${species_reference_sa_path} .
    """
}

process GET_TARGET_INFORMATION {
    cache 'lenient'
    publishDir "${params.outdir}/full_probe_info/${sample_name}/"
    input:
        tuple val(sample_name),val(species),path(R1), path(R2), val(attb_name), val(attp_name),val(umi_type),val(probes_name),val(cargo),val(group),path(reference_genome), path(reference_index1),path(reference_index2),path(reference_index3),path(reference_index4),path(reference_index5),path(reference_index6)
        path cosmic_info
        path gtex_info
        val benchling_warehouse_username
        val benchling_warehouse_password
        val benchling_warehouse_url
        val benchling_sdk_api_key
        val benchling_api_url

    output:
        tuple val(sample_name), val(group), path("${sample_name}_target_info.csv"), emit: target_info
        tuple val(sample_name), val(group), val(cargo), path("cargo.fasta"), emit: cargo_reference

    script:
        """
        export WAREHOUSE_USERNAME='${benchling_warehouse_username}'
        export WAREHOUSE_PASSWORD='${benchling_warehouse_password}'
        export WAREHOUSE_URL='${benchling_warehouse_url}'
        export API_KEY='${benchling_sdk_api_key}'
        export API_URL='${benchling_api_url}'
        get_target_info.py --probes_name "${probes_name}" --cosmic_info "${cosmic_info}" --gtex_info "${gtex_info}" --attp_name "${attp_name}" --reference "${reference_genome}" --cargo "${cargo}" --sample_name "${sample_name}"
        """
}

process GENERATE_REFERENCE_CARGO_GENOME {
    input:
        tuple val(cargo_name), path(reference_fasta)

    output:
        tuple val(cargo_name), path("*.{fa,fna,fasta}"), emit: reference_fasta
        tuple val(cargo_name), path("*.{fa,fna,fasta}.*"), emit: reference_index

    script:
    """
    python -c 'import sys; sys.path.append("${workflow.projectDir}/bin/"); import get_target_info;  get_target_info.download_cargo_genome("${cargo_name}")'
    cat ${reference_fasta} cargo.fasta > ref_${cargo_name}.fa
    bwa index -@ 32 ref_${cargo_name}.fa
    samtools faidx ref_${cargo_name}.fa
    """
}

process ADAPTER_AND_POLY_G_TRIM {
    cache 'lenient'
    publishDir "${params.outdir}/raw_fastq/${sample_name}/", pattern: '*.fastq.gz'
    publishDir "${params.outdir}/trimmed_fastq/${sample_name}/", pattern: '*trimmed*.fastq.gz'
    publishDir "${params.outdir}/input_probe_sheet/${sample_name}/", pattern: '*.csv'

    input:
        tuple val(sample_name),val(species),path(R1), path(R2), val(attb_name), val(attp_name),val(umi_type),val(probes_name),val(cargo),val(group)
    output:
        path "${R1}"
        path "${R2}"
        tuple val(sample_name), val(group), val(umi_type), path("${sample_name}_trimmed.fastq.gz"), emit: trimmed_fastq
        tuple val(sample_name), val(group), path("${sample_name}_fastp.json"), emit: fastp_stats
    script:
        umi_loc = ''
        umi_len = ''
        umi_skip = ''
        umi_type = "Twist"

        if (umi_type == "Twist") {
            umi_loc = 'per_read'
            umi_len = '5'
            umi_skip = '2'
            umi_params = "--umi_loc=${umi_loc} --umi_len=${umi_len} --umi_skip=${umi_skip}"
        } else if (umi_type == "LMPCR" || umi_type == "LM-PCR") {
            umi_loc = 'read1'
            umi_len = '11'
            umi_params = "--umi_loc=${umi_loc} --umi_len=${umi_len}"
        } else if (umi_type == "xGen") {
            umi_loc = 'read1'
            umi_len = '9'
            umi_params = "--umi_loc=${umi_loc} --umi_len=${umi_len}"
        } else if (umi_type == "None") {
            umi_params = ""
        }
        """
        fastp -m -c --dont_eval_duplication --disable_adapter_trimming --low_complexity_filter --overlap_len_require 10 -i ${R1} -I ${R2} --merged_out ${sample_name}_trimmed.fastq.gz --include_unmerged -w 16 -g -j ${sample_name}_fastp.json -U ${umi_params} --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
        """
}

process ALIGN_READS {
    cache 'lenient'
    publishDir "${params.outdir}/initial_alignments/${sample_name}/", pattern: '*_initial_alignment.bam*'
    publishDir "${params.outdir}/deduped_alignments/${sample_name}/", pattern: '*_deduped_alignment.bam*'

    input:
        tuple val(sample_name), val(group), val(umi_type), path(fastq), val(cargo), path(genome_reference), path(reference_index)
        val initial_mapper
    output:
        tuple val(sample_name), val(group), path("${sample_name}_initial_alignment.bam"), path("${sample_name}_initial_alignment.bam.bai"), emit: original_alignment_bam
        tuple val(sample_name), val(group), path("${sample_name}_deduped_alignment.bam"), path("${sample_name}_deduped_alignment.bam.bai"), emit: deduped_alignment_bam

    script:

    def alignment_command = ""
    if (initial_mapper == "minimap2") {
        alignment_command = "minimap2 -ax sr -t 96 ${genome_reference} ${fastq}"
    } else if (initial_mapper == "bwa") {
        alignment_command = "bwa mem -p -t 32 ${genome_reference} ${fastq}"
    } else {
        error "Unsupported initial_mapper: $initial_mapper"
    }
    if (umi_type != "None") {
        """
        # align in the single-end & paired-end mix mode -p
        $alignment_command > ${sample_name}_initial_alignment.sam
        samtools sort ${sample_name}_initial_alignment.sam -@ 32 -o ${sample_name}_initial_alignment.bam
        samtools index ${sample_name}_initial_alignment.bam

        # dedup in single-end mode
        umi_tools dedup -I ${sample_name}_initial_alignment.bam --umi-separator ":" -S ${sample_name}_deduped_alignment.bam --method unique
        samtools index ${sample_name}_deduped_alignment.bam

        rm *.sam
        """
    } else {
        """
        $alignment_command | sambamba view -S -t 96 -f bam -o ${sample_name}_initial_alignment.bam /dev/stdin
        sambamba sort --show-progress -t 96 -o ${sample_name}_initial_alignment_sorted.bam ${sample_name}_initial_alignment.bam
        mv ${sample_name}_initial_alignment_sorted.bam ${sample_name}_initial_alignment.bam
        samtools index ${sample_name}_initial_alignment.bam
        cp ${sample_name}_initial_alignment.bam ${sample_name}_deduped_alignment.bam
        rm ${sample_name}_initial_alignment.sam
        samtools index ${sample_name}_deduped_alignment.bam
        """
    }
}

process EXTRACT_TARGET_READS {
    cache 'lenient'
    publishDir "${params.outdir}/probe_specific_reads/${sample_name}/"
    input:
        tuple val(sample_name), val(group), path(target_info), path(bam_file), path(bam_file_index)
    output:
        tuple val(sample_name), val(group), path("*.fastq"), emit: extracted_reads
        tuple val(sample_name), val(group), path("${sample_name}_read_counts_per_site.csv"), emit: read_counts_per_site_file

    script:
    """
    extract_target_reads.py --target_info ${target_info} --bam_file ${bam_file} --sample_name ${sample_name}
    """
}

process GENERATE_AMPLICONS {
    cache 'lenient'
    publishDir "${params.outdir}/generated_amplicons/${sample_name}/"
    input:
        tuple val(sample_name), val(group), path(target_info)
    output:
        tuple val(sample_name), val(group), path("*.fasta")
    script:
    """
    create_amplicon_files.py --target_info ${target_info}
    """
}

process ALIGN_TARGET_READS {
    cache 'lenient'
    publishDir "${params.outdir}/read_to_probe_alignment/${sample_name}/", pattern:'*.bam*'
    input:
        tuple val(sample_name), val(group), path(target_info), path(fastq_files), path(amplicon_files)
    output:
        tuple val(sample_name), val(group), path("*.bam*")
    script:
    """
    align_extracted_reads.py --target_info ${target_info} --fastq_files ${fastq_files} --amplicon_files ${amplicon_files} --sample_name ${sample_name}
    """
}

process MEASURE_INTEGRATION {
    cache 'lenient'
    publishDir "${params.outdir}/integration_stats_tables/${sample_name}/", pattern: '*integration_stats.csv'
    input:
        tuple val(sample_name), val(group), path(target_info), path(alignments)
    output:
        tuple val(sample_name), val(group), path("${sample_name}_integration_stats.csv"), emit: integration_stats_file
        tuple val(sample_name), val(group), path("*.fastq"), emit: edited_reads, optional: true
    script:
    """
    compute_integration_percentage.py --target_info ${target_info} --bam ${alignments} --sample_name ${sample_name}
    """
}

process UPDATE_BENCHLING_WITH_VALIDATED_SITES {
    input:
        val project_id
        tuple val(sample_name), val(group), path(integration_csv)
    output:

    script:
    """
    update_benchling_with_offtargets.py --project_id ${project_id} --integration_csv ${integration_csv}
    """
}

process GATHER_QC_INFO {
    cache 'lenient'
    publishDir "${params.outdir}/qc_info/${sample_name}/"
    input:
        tuple val(sample_name), val(group), path(json_file),path(original_bam_file), path(original_bam_file_index),path(deduped_bam_file), path(deduped_bam_file_index)
    output:
        tuple val(sample_name), val(group), path("${sample_name}_qc_summary.csv"), emit: qc_summary_file
    script:
    """
    gather_qc_stats.py --json_file ${json_file} --sample_name ${sample_name} --original_bam_file ${original_bam_file} --deduped_bam_file ${deduped_bam_file}
    """
}

process ALIGNMENT_VISUALIZATION {
    cache 'lenient'
    publishDir "${params.outdir}/alignment_visualizations/${sample_name}/attL_alignments/", pattern:'*attL*.html'
    publishDir "${params.outdir}/alignment_visualizations/${sample_name}/attR_alignments/", pattern:'*attR*.html'
    publishDir "${params.outdir}/alignment_visualizations/${sample_name}/beacon_alignments/", pattern:'*beacon*.html'
    publishDir "${params.outdir}/alignment_visualizations/${sample_name}/wt_alignments/", pattern:'*wt*.html'

    input:
        tuple val(sample_name), val(group), path(target_info), path(amplicons), path(probe_read_alignments), path(edited_reads)
        val(bam2html_path)
    output:
        val(sample_name), emit: sample_name
        path("*.html"), optional: true

    script:
    """
    alignment_visualization.py --bam2html_path ${bam2html_path} --bam ${probe_read_alignments} --target_info ${target_info} --sample_name ${sample_name}
    """
}

process GENERATE_REPORT {
    cache 'lenient'
    publishDir "${params.outdir}"
    input:
        path project_config_file
        path integration_stats_files
        path read_counts_per_site_files
        path qc_summary_files
        path extracted_reads
        val collapse_condition
        val project_name
    output:
        path "*.xlsx", emit: excel_output

    script:
        """
        ulimit -s 65536
        collate_results.py --project_config_file ${project_config_file} --integration_stats_files ${integration_stats_files} --read_counts_per_site_files ${read_counts_per_site_files} --qc_summary_files ${qc_summary_files} --extracted_reads_dirs ${extracted_reads} --collapse_condition ${collapse_condition} --project_name ${project_name}
        """
    }

process MULTIQC {
    cache 'lenient'
    publishDir "${params.outdir}"
    input:
        path fastp_jsons
    output:
        file "multiqc_report.html"
        file "multiqc_data"

    script:
    """
        multiqc .
    """
}

process CREATE_PYTHON_NOTEBOOK_REPORT {
    cache 'lenient'
    publishDir "${params.outdir}"

    input:
        path excel_file
        val notebook_template
    output:
        path 'report.html'
    script:
    """
    papermill ${notebook_template} report.ipynb -p results_file ${excel_file}
    jupyter nbconvert --to html --no-input report.ipynb
    """
}

process CREATE_QUILT_PACKAGE {
    input:
        val output_folder
        path notebook_report
        path cas_bed
        val project_id
        val bucket_name
        val quilt_output
        val benchling_warehouse_username
        val benchling_warehouse_password
        val benchling_warehouse_url
        val benchling_sdk_api_key
        val benchling_api_url
    output:

    script:
    """
    export WAREHOUSE_USERNAME='${benchling_warehouse_username}'
    export WAREHOUSE_PASSWORD='${benchling_warehouse_password}'
    export WAREHOUSE_URL='${benchling_warehouse_url}'
    export API_KEY='${benchling_sdk_api_key}'
    export API_URL='${benchling_api_url}'
    create_quilt_package.py --output_folder ${workflow.launchDir}/${output_folder} --project_id ${project_id} --bucket_name ${bucket_name} --package_name ${quilt_output}
    """
}

process TRANSLOCATION_DETECTION {
    cache 'lenient'
    publishDir "${params.outdir}/translocation/"

    input:
        tuple val(sample_name), val(group), path(target_info), val(cargo_name), path(reference_genome), path(reference_index), path(bam_file), path(bam_file_index)

    output:
        tuple val(sample_name), val(group), path("*.bnd.bed"), path("*.cargo.bed"), emit: bnd
        path '*.vcf*'

    script:
    """
    # delly SV call
    delly call -g ${reference_genome} ${bam_file} -o ${sample_name}.delly.vcf -q 0 -r 0 -c 1 -z 0 -m 0 -t DEL,INV,BND &> log
    bcftools query ${sample_name}.delly.vcf -i 'SVTYPE="BND"' -f "%CHROM\\t%POS\\t%CHR2\\t%POS2\\t%ID\\t%PE\\t%SR\\n" > ${sample_name}.bnd.bed
    bcftools query ${sample_name}.delly.vcf -i 'SVTYPE="INV" || SVTYPE="DEL"' -f "%CHROM\\t%POS\\t%CHROM\\t%END\\t%ID\\t%PE\\t%SR\\n" >> ${sample_name}.bnd.bed

    # extract the paired-end reads, one end mapped to reference, and one end to the cargo
    samtools view -f1 -F2304 ${bam_file} -e '(rname=="$cargo_name" || rnext=="$cargo_name") && rname!=rnext' -b -@ 10 > ${sample_name}.initial_alignment.cargo.bam
    samtools sort -n ${sample_name}.initial_alignment.cargo.bam | bedtools bamtobed -bedpe | awk -F "\\t" '{ OFS="\\t"; split(\$7,f,":"); print \$4,\$5,\$6,f[length(f)] }' | sort -u | sort -k1,1 -k2,2n > ${sample_name}.cargo.bed
    """
}

process INTERSECT_CAS_DATABASE {
    publishDir "${params.outdir}/translocation/"

    input:
        val dinucleotides
        tuple val(sample_name), val(group), path(bnd_file), path(cargo_file), val(species)
    output:
        path '*.cas.bed', emit: cas_bed
        path '*.cargo.cas.csv', emit: cargo_cas_bed

    script:
    def args = dinucleotides == "" ? "": "--dinucleotides ${dinucleotides}"
    """
    get_cas_info.py --species $species $args | sort -k1,1 -k2,2n | bedtools groupby -g 1,2,3 -c 4,5,6 -o distinct,distinct,distinct > CAS.cut.bed
    awk -F "\\t" '{OFS="\\t"; print \$1,\$2-1,\$2,\$5,\$6+\$7"\\n"\$3,\$4-1,\$4,\$5,\$6+\$7}' ${bnd_file} | sort -k1,1 -k2,2n | bedtools closest -a stdin -b CAS.cut.bed -d | awk -F "\\t" '{OFS="\\t"; if(\$12<=10 && \$12>=0) print}' > ${sample_name}.bnd.cas.bed

    bedtools closest -a ${cargo_file} -b CAS.cut.bed -d | awk -F "\\t" '{OFS="\\t"; if(\$11>0 && \$11<=100) print }' > tmp
    if [ -s tmp ]
    then
        sort -k8 tmp | bedtools groupby -g 8 -c 4 -o count > ${sample_name}.cargo.cas.csv
    else
        touch ${sample_name}.cargo.cas.csv
    fi
    """
}

workflow {

    log.info """\
        #################################################
        #                                               #
        #            Hybrid Capture Analysis (HCA)      #
        #                                               #
        #                                               #
        #################################################
         ${params.manifest.name} v${params.manifest.version}
         ${params.manifest.author}
         ${params.manifest.description}
         ==========================
         project name : ${params.project_id}
         input from   : ${params.samplesheet}
         output to    : ${params.outdir}

         other parameters:
         reference genome: ${params.reference}
         reference genome index: ${params.reference_index_location}
         initial mapper: ${params.initial_mapper}
         attP, left side: ${params.ATTP_REG}
         attP, right side: ${params.ATTP_PRIME}
         collapse condition: ${params.collapse_condition}
         --
         run as       : ${workflow.commandLine}
         started at   : ${workflow.start}
         config files : ${workflow.configFiles}
         container    : ${workflow.containerEngine}:${workflow.container}
         """
         .stripIndent()

    raw_reads = DOWNLOAD_READS(params.project_id)
    samplesheet = GET_PROJECT_INFO(params.project_id,raw_reads,params.BENCHLING_WAREHOUSE_USERNAME,params.BENCHLING_WAREHOUSE_PASSWORD,params.BENCHLING_WAREHOUSE_URL,params.BENCHLING_API_KEY,params.BENCHLING_API_URL)


    // Existing logic to handle FASTQ file input
    samplesheet
    .splitCsv(header: true, sep: ',')
    .map { row ->
        tuple(
            row.sample_name,
            row.species,
            file(row.read1),
            file(row.read2),
            row.attb,
            row.attp,
            row.umi_type,
            row.probes_name,
            row.cargo_name,
            row.group
        )
    }
    .set { input_ch }

    input_ch.map { tuple -> tuple[1] }  // Select the second element from each tuple
        .unique()                   // Remove duplicates to get unique items
        .set {unique_species_ch}                     // Print the unique items

    reference_genome = DOWNLOAD_REFERENCE_GENOME(unique_species_ch)
    gtex_data = DOWNLOAD_GTEX_DATA()
    trimmed_and_merged_fastq = ADAPTER_AND_POLY_G_TRIM(input_ch)

    input_ch
        .combine(reference_genome.reference_fasta)
        .combine(reference_genome.reference_index)
        .set{probe_info_input_ch}

    probe_information = GET_TARGET_INFORMATION(probe_info_input_ch,params.cosmic_info,gtex_data,params.BENCHLING_WAREHOUSE_USERNAME,params.BENCHLING_WAREHOUSE_PASSWORD,params.BENCHLING_WAREHOUSE_URL,params.BENCHLING_API_KEY,params.BENCHLING_API_URL)

    amplicon_files = GENERATE_AMPLICONS(probe_information.target_info)

    // generate ref and cargo reference
    probe_information.cargo_reference.map { it[2] }
        .unique()
        .combine(reference_genome.reference_fasta)
        .set { cargo_ch }
    cargo_ch.view()
    reference_cargo_genome = GENERATE_REFERENCE_CARGO_GENOME(cargo_ch)

    probe_information.cargo_reference
        .map { [it[2], it[0], it[1]] }
        .combine(reference_cargo_genome.reference_fasta, by:[0])
        .combine(reference_cargo_genome.reference_index, by:[0])
        .map { [it[1], it[2], it[0], it[3], it[4]]}
        .set { ref_ch }

    trimmed_and_merged_fastq.trimmed_fastq
        .combine(ref_ch, by:[0,1])
        .set{align_reads_input_ch}
    initial_alignment = ALIGN_READS(align_reads_input_ch, params.initial_mapper)

    probe_information.target_info
        .combine(initial_alignment.deduped_alignment_bam, by: [0,1])
        .set{target_info_and_deduped_alignment_ch}

    extract_target_reads_out = EXTRACT_TARGET_READS(target_info_and_deduped_alignment_ch)

    probe_information.target_info
        .combine(extract_target_reads_out.extracted_reads, by: [0,1])
        .combine(amplicon_files, by:[0,1])
        .set{align_target_reads_input_ch}

    align_target_reads_out = ALIGN_TARGET_READS(align_target_reads_input_ch)

    probe_information.target_info
        .combine(align_target_reads_out,by:[0,1])
        .set{measure_integration_input_ch}

    measure_integration_out = MEASURE_INTEGRATION(measure_integration_input_ch)

    // UPDATE_BENCHLING_WITH_VALIDATED_SITES(params.project_id,measure_integration_out.integration_stats_file)

    trimmed_and_merged_fastq.fastp_stats
        .combine(initial_alignment.original_alignment_bam,by:[0,1])
        .combine(initial_alignment.deduped_alignment_bam,by:[0,1])
        .set{gather_qc_info_input_ch}
    qc_summary = GATHER_QC_INFO(gather_qc_info_input_ch)

    probe_information.target_info
        .combine(amplicon_files,by:[0,1])
        .combine(align_target_reads_out,by:[0,1])
        .combine(measure_integration_out.edited_reads,by:[0,1])
        .set{alignment_viz_input_ch}
    ALIGNMENT_VISUALIZATION(alignment_viz_input_ch, params.bam2html_path)

    measure_integration_out.integration_stats_file
        .collect(flat:false)
        .flatMap{ it }
        .map{ tuple -> tuple[2] }
        .collect()
        .set{integration_stats_files_ch}

    extract_target_reads_out.read_counts_per_site_file
        .collect(flat:false)
        .flatMap{ it }
        .map{ tuple -> tuple[2]}
        .collect()
        .set{read_counts_per_site_files_ch}

    qc_summary.qc_summary_file
        .collect(flat:false)
        .flatMap{ it }
        .map{ tuple -> tuple[2]}
        .collect()
        .set{qc_summary_files_ch}

    extract_target_reads_out.extracted_reads
        .collect(flat:false)
        .flatMap{ it }
        .map{ tuple -> tuple[2]}
        .collect()
        .set{extracted_reads_files_ch}

    report_excel_file = GENERATE_REPORT(samplesheet,integration_stats_files_ch,read_counts_per_site_files_ch,qc_summary_files_ch,extracted_reads_files_ch, params.collapse_condition, params.project_id)

    // ** MULTIQC REPORT **
    trimmed_and_merged_fastq.fastp_stats
        .flatten()  // Flatten the list
        .filter { it.toString().endsWith('.json') }  // Filter out only the paths ending with .json
        .collect()
        .ifEmpty([])
        .set{multiqc_input_ch}
    MULTIQC(multiqc_input_ch)

    // ** TRANSLOCATION DETECTION **
    probe_information.target_info
        .combine(ref_ch, by: [0,1])
        .combine(initial_alignment.original_alignment_bam, by: [0,1])
        .set{translocation_detection_input_ch}
    TRANSLOCATION_DETECTION(translocation_detection_input_ch)

    TRANSLOCATION_DETECTION.out.bnd
        .combine(unique_species_ch)
        .set { bnd_ch }
    bnd_ch.view()
    intersect_cas_database_out = INTERSECT_CAS_DATABASE(params.dinucleotides, bnd_ch)

    // ** CREATE HTML REPORT **

    html_report = CREATE_PYTHON_NOTEBOOK_REPORT(report_excel_file, params.notebook_template)

    // CREATE_QUILT_PACKAGE(params.outdir,html_report,intersect_cas_database_out.cas_bed.collect(),params.project_id,params.bucket_name,params.quilt_package_name,params.BENCHLING_WAREHOUSE_USERNAME,params.BENCHLING_WAREHOUSE_PASSWORD,params.BENCHLING_WAREHOUSE_URL,params.BENCHLING_API_KEY,params.BENCHLING_API_URL)

}
