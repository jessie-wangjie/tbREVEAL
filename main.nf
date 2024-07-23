nextflow.enable.dsl=2


params.outdir = "output"
params.initial_mapper = 'bwa'
params.notebook_template = "${workflow.projectDir}/bin/report_generation.ipynb"
params.bam2html_path = "${workflow.projectDir}/bin/utils/bam2html.py"
params.dinucleotides = ''
params.project_id=''


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

process ADAPTER_AND_POLY_G_TRIM {
    cache 'lenient'
    publishDir "${params.outdir}/raw_fastq/${sample_name}/", pattern: '*.fastq.gz', overwrite: true
    publishDir "${params.outdir}/trimmed_fastq/${sample_name}/", pattern: '*trimmed*.fastq.gz', overwrite: true
    publishDir "${params.outdir}/input_probe_sheet/${sample_name}/", pattern: '*.csv', overwrite: true

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
    publishDir "${params.outdir}/initial_alignments/${sample_name}/", pattern: '*_initial_alignment.bam*', overwrite: true
    publishDir "${params.outdir}/deduped_alignments/${sample_name}/", pattern: '*_deduped_alignment.bam*', overwrite: true

    input:
        tuple val(sample_name), val(group), val(umi_type), path(fastq), path(genome_reference), path(reference_index1),path(reference_index2),path(reference_index3),path(reference_index4),path(reference_index5),path(reference_index6),path(cargo_reference)
        val initial_mapper
    output:
        tuple val(sample_name), val(group), path("${sample_name}_initial_alignment.bam"), path("${sample_name}_initial_alignment.bam.bai"), emit: original_alignment_bam
        tuple val(sample_name), val(group), path("${sample_name}_deduped_alignment.bam"), path("${sample_name}_deduped_alignment.bam.bai"), emit: deduped_alignment_bam

    script:

    def alignment_command = ""
    if (initial_mapper == "minimap2") {
        alignment_command = "minimap2 -ax sr -t 96 ${genome_reference} ${fastq}"
    } else if (initial_mapper == "bwa") {
        alignment_command = "bwa mem -t 96 ${genome_reference} ${fastq}"
        cargo_alignment_command = "bwa mem -t 96 ${cargo_reference}"
    } else {
        error "Unsupported initial_mapper: $initial_mapper"
    }
    if (umi_type != "None") {
        """
        $alignment_command | sambamba view -S -t 96 -f bam -o ${sample_name}_initial_alignment.bam /dev/stdin
        sambamba sort --show-progress -t 96 -o ${sample_name}_initial_alignment_sorted.bam ${sample_name}_initial_alignment.bam
        mv ${sample_name}_initial_alignment_sorted.bam ${sample_name}_initial_alignment.bam
        mv ${sample_name}_initial_alignment_sorted.bam.bai ${sample_name}_initial_alignment.bam.bai
        samtools index ${sample_name}_initial_alignment.bam
        umi_tools dedup -I ${sample_name}_initial_alignment.bam --umi-separator ":" -S ${sample_name}_deduped_alignment.bam --method unique
        samtools index ${sample_name}_deduped_alignment.bam

        samtools view -f 4 ${sample_name}_initial_alignment.bam | samtools bam2fq | pigz -p 8 > ${sample_name}_unaligned_reads.fastq.gz
        bwa index ${cargo_reference}
        $cargo_alignment_command ${sample_name}_unaligned_reads.fastq.gz | sambamba view -S -t 96 -f bam -o ${sample_name}_unaligned_reads_aligned_to_cargo.bam /dev/stdin
        sambamba sort --show-progress -t 96 -o ${sample_name}_unaligned_reads_aligned_to_cargo_sorted.bam ${sample_name}_unaligned_reads_aligned_to_cargo.bam
        mv ${sample_name}_unaligned_reads_aligned_to_cargo_sorted.bam ${sample_name}_unaligned_reads_aligned_to_cargo.bam
        mv ${sample_name}_unaligned_reads_aligned_to_cargo_sorted.bam.bai ${sample_name}_unaligned_reads_aligned_to_cargo.bam.bai
        samtools index ${sample_name}_unaligned_reads_aligned_to_cargo.bam
        umi_tools dedup -I ${sample_name}_unaligned_reads_aligned_to_cargo.bam --umi-separator ":" -S ${sample_name}_unaligned_reads_deduped_aligned_to_cargo.bam --method unique
        samtools index ${sample_name}_unaligned_reads_deduped_aligned_to_cargo.bam

        samtools merge -o final.bam ${sample_name}_deduped_alignment.bam ${sample_name}_unaligned_reads_deduped_aligned_to_cargo.bam
        mv final.bam ${sample_name}_deduped_alignment.bam
        samtools index ${sample_name}_deduped_alignment.bam
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
    publishDir "${params.outdir}/probe_specific_reads/${sample_name}/", overwrite: true
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
    publishDir "${params.outdir}/generated_amplicons/${sample_name}/", overwrite: true
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
    publishDir "${params.outdir}/read_to_probe_alignment/${sample_name}/", pattern:'*.bam*', overwrite: true
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
    publishDir "${params.outdir}/integration_stats_tables/${sample_name}/", pattern: '*integration_stats.csv', overwrite: true
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

process GATHER_QC_INFO {
    cache 'lenient'
    publishDir "${params.outdir}/qc_info/${sample_name}/", overwrite: true
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
    publishDir "${params.outdir}/alignment_visualizations/${sample_name}/", pattern:'*.html', overwrite: true

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
    publishDir "${params.outdir}", overwrite: true
    input:
        path project_config_file
        path integration_stats_files
        path read_counts_per_site_files
        path qc_summary_files
        path extracted_reads
        val project_name
    output:
        path "*.xlsx", emit: excel_output

    script:
        """
        ulimit -s 65536
        collate_results.py --project_config_file ${project_config_file} --integration_stats_files ${integration_stats_files} --read_counts_per_site_files ${read_counts_per_site_files} --qc_summary_files ${qc_summary_files} --extracted_reads_dirs ${extracted_reads} --project_name "${project_name}"
        """
    }

process MULTIQC {
    cache 'lenient'
    publishDir "${params.outdir}", overwrite: true
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
    publishDir "${params.outdir}", overwrite: true

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


process TRANSLOCATION_DETECTION {
    cache 'lenient'
    publishDir "${params.outdir}/translocation/", overwrite: true
    input:
        tuple val(sample_name), val(group),path(target_info), path(cargo), path(bam_file), path(bam_file_index), path(reference_genome)
    output:
        tuple val(sample_name), val(group), path("*.bnd.bed"), emit: bnd
        path '*.vcf*'
    script:
    """
    cat ${reference_genome} ${cargo} > combined_genomes.fasta
    samtools faidx combined_genomes.fasta
    delly call -g combined_genomes.fasta ${bam_file} -o ${sample_name}.delly.vcf -q 0 -r 0 -c 1 -z 0 -m 0 -t DEL,INV,BND &> log
    bcftools query ${sample_name}.delly.vcf -i 'SVTYPE=\"BND\"' -f "%CHROM\\t%POS\\t%CHR2\\t%POS2\\t%ID\\t%PE\\t%SR\\n" > ${sample_name}.bnd.bed
    bcftools query ${sample_name}.delly.vcf -i 'SVTYPE=\"INV\" || SVTYPE=\"DEL\"' -f "%CHROM\\t%POS\\t%CHROM\\t%END\\t%ID\\t%PE\\t%SR\\n" >> ${sample_name}.bnd.bed
    """
}

process INTERSECT_CAS_DATABASE {
    publishDir "${params.outdir}/translocation/", overwrite: true

    input:
        val dinucleotides
        tuple val(sample_name), val(group), path(bnd_file)
        val benchling_warehouse_username
        val benchling_warehouse_password
        val benchling_warehouse_url
        val benchling_sdk_api_key
        val benchling_api_url
    output:
        path '*.cas.bed', emit: cas_bed

    script:
    def args = dinucleotides == "" ? "": "--dinucleotides ${dinucleotides}"
    """
    export WAREHOUSE_USERNAME='${benchling_warehouse_username}'
    export WAREHOUSE_PASSWORD='${benchling_warehouse_password}'
    export WAREHOUSE_URL='${benchling_warehouse_url}'
    export API_KEY='${benchling_sdk_api_key}'
    export API_URL='${benchling_api_url}'
    get_cas_info.py $args | sort -k1,1 -k2,2n | bedtools groupby -g 1,2,3 -c 4,5,6 -o distinct,distinct,distinct > CAS.cut.bed
    awk -F \"\\t\" '{OFS=\"\\t\"; print \$1,\$2-1,\$2,\$5,\$6+\$7\"\\n\"\$3,\$4-1,\$4,\$5,\$6+\$7}' ${bnd_file} | sort -k1,1 -k2,2n | bedtools closest -a stdin -b CAS.cut.bed -d | awk -F \"\\t\" '{OFS=\"\\t\"; if(\$12<=10 && \$12>=0) print}' > ${sample_name}.bnd.cas.bed
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
         --
         run as       : ${workflow.commandLine}
         started at   : ${workflow.start}
         config files : ${workflow.configFiles}
         container    : ${workflow.containerEngine}:${workflow.container}
         """
         .stripIndent()

    raw_reads = DOWNLOAD_READS(params.project_id,params.BS_ACCESS_TOKEN,params.BS_API_SERVER)
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

    reference_genome = DOWNLOAD_REFERENCE_GENOME(unique_species_ch,params.AWS_ACCESS_KEY_ID,params.AWS_SECRET_ACCESS_KEY)
    gtex_data = DOWNLOAD_GTEX_DATA()
    trimmed_and_merged_fastq = ADAPTER_AND_POLY_G_TRIM(input_ch)

    input_ch
        .combine(reference_genome.reference_fasta)
        .combine(reference_genome.reference_index)
        .set{probe_info_input_ch}

    probe_information = GET_TARGET_INFORMATION(probe_info_input_ch,params.cosmic_info,gtex_data,params.BENCHLING_WAREHOUSE_USERNAME,params.BENCHLING_WAREHOUSE_PASSWORD,params.BENCHLING_WAREHOUSE_URL,params.BENCHLING_API_KEY,params.BENCHLING_API_URL)

    amplicon_files = GENERATE_AMPLICONS(probe_information.target_info)

    trimmed_and_merged_fastq.trimmed_fastq
        .combine(reference_genome.reference_fasta)
        .combine(reference_genome.reference_index)
        .combine(probe_information.cargo_reference,by:[0,1])
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

    UPDATE_BENCHLING_WITH_VALIDATED_SITES(params.project_id,measure_integration_out.integration_stats_file,params.BENCHLING_WAREHOUSE_USERNAME,params.BENCHLING_WAREHOUSE_PASSWORD,params.BENCHLING_WAREHOUSE_URL,params.BENCHLING_API_KEY,params.BENCHLING_API_URL)

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

    integration_stats_files_ch.view()

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

    report_excel_file = GENERATE_REPORT(samplesheet,integration_stats_files_ch,read_counts_per_site_files_ch,qc_summary_files_ch,extracted_reads_files_ch, params.project_id)

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
        .combine(probe_information.cargo_reference, by: [0,1])
        .combine(initial_alignment.deduped_alignment_bam, by: [0,1])
        .combine(reference_genome.reference_fasta)
        .set{translocation_detection_input_ch}

    TRANSLOCATION_DETECTION(translocation_detection_input_ch)

    intersect_cas_database_out = INTERSECT_CAS_DATABASE(params.dinucleotides, TRANSLOCATION_DETECTION.out.bnd,params.BENCHLING_WAREHOUSE_USERNAME,params.BENCHLING_WAREHOUSE_PASSWORD,params.BENCHLING_WAREHOUSE_URL,params.BENCHLING_API_KEY,params.BENCHLING_API_URL)

    // ** CREATE HTML REPORT **

    html_report = CREATE_PYTHON_NOTEBOOK_REPORT(report_excel_file, params.notebook_template)

    CREATE_QUILT_PACKAGE(params.outdir,html_report,intersect_cas_database_out.cas_bed.collect(),params.project_id,params.bucket_name,params.quilt_package_name,params.BENCHLING_WAREHOUSE_USERNAME,params.BENCHLING_WAREHOUSE_PASSWORD,params.BENCHLING_WAREHOUSE_URL,params.BENCHLING_API_KEY,params.BENCHLING_API_URL,params.AWS_ACCESS_KEY_ID,params.AWS_SECRET_ACCESS_KEY)

}
