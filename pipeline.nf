nextflow.enable.dsl=2

params.reference = "/data/references/hg38.fa"
params.outdir = ""
params.ATTP_REG = 'GTGGTTTGTCTGGTCAACCACCGCGGT'
params.ATTP_PRIME = 'CTCAGTGGTGTACGGTACAAACCCA'
params.deduplication_method = 'STRICT'
params.collapse_condition = 'Complete'
params.project_name = ''
params.samplesheet = ''

process ADAPTER_AND_POLY_G_TRIM {
    cache 'lenient'
    maxForks 1
    publishDir "${params.outdir}/${sample_name}/qc/", pattern: '*fastp.json'
    publishDir "${params.outdir}/${sample_name}/raw_fastq/", pattern: '*.fastq.gz'
    publishDir "${params.outdir}/${sample_name}/trimmed_fastq/", pattern: '*.fastq.bgz'
    publishDir "${params.outdir}/${sample_name}/probe_information_table/", pattern: '*.csv'

    input:
        tuple val(sample_name), path(R1), path(R2), path(metadata)
    output:
        path "${R1}"
        path "${R2}"
        path "${metadata}"
        tuple val(sample_name), path("${sample_name}_trimmed.fastq.bgz"), emit: trimmed_fastq
        tuple val(sample_name), path("${sample_name}_fastp.json"), emit: fastp_stats

    script:
        """
        fastp -m -c --include_unmerged --low_complexity_filter --cut_right --overlap_len_require 10 -i ${R1} -I ${R2} --merged_out ${sample_name}_trimmed.fastq.gz -w 16 -g -j ${sample_name}_fastp.json
        gunzip -c ${sample_name}_trimmed.fastq.gz | bgzip -@ 8 > ${sample_name}_trimmed.fastq.bgz
        """
}

process CREATE_MINIMAP2_REFERENCE_INDEX {
    cache 'lenient'

    input:
        path reference
    output:
        path 'target.mmi'

    script:
        """
        minimap2 -x sr -d target.mmi ${reference}
        """
}

process ALIGN_READS {
    cache 'lenient'
    maxForks 1
    // memory 1.GB
    publishDir "${params.outdir}/${sample_name}/"

    input:
        //tuple val(sample_name), path(R1), path(R2), path(metadata)
        path reference
        path reference_index
        tuple val(sample_name), path(fastq)
        val dedup_method

    output:
        val(sample_name), emit: sample_name
        tuple path("${sample_name}_initial_alignment.bam"), path("${sample_name}_initial_alignment.bam.bai"), emit: initial_alignment_bam
        path "${reference}.fai"

    script:
    if (dedup_method == "STRICT")
        """
        minimap2 -ax sr -t 16 ${reference_index} ${fastq} > ${sample_name}_initial_alignment.sam 
        samtools view -@ 16 -b ${sample_name}_initial_alignment.sam > ${sample_name}_initial_alignment.bam
        samtools sort -@ 16 ${sample_name}_initial_alignment.bam > ${sample_name}_initial_alignment_sorted.bam
        mv ${sample_name}_initial_alignment_sorted.bam ${sample_name}_initial_alignment.bam
        rm ${sample_name}_initial_alignment.sam
        samtools sort -@ 16 -n -o ${sample_name}_namesort.bam ${sample_name}_initial_alignment.bam
        samtools fixmate -m ${sample_name}_namesort.bam ${sample_name}_fixmate.bam
        samtools sort -@ 16  -o ${sample_name}_positionsort.bam ${sample_name}_fixmate.bam
        samtools markdup --barcode-rgx ":(\\w+)\$"  -r ${sample_name}_positionsort.bam ${sample_name}_markdup.bam  
        mv ${sample_name}_markdup.bam ${sample_name}_initial_alignment.bam
        samtools index ${sample_name}_initial_alignment.bam
        samtools faidx ${reference}
        rm ${sample_name}_positionsort.bam
        rm ${sample_name}_fixmate.bam
        rm ${sample_name}_namesort.bam
        """
    else if (dedup_method == "LOOSE")
        """
        minimap2 -ax sr -t 16 ${reference_index} ${fastq} > ${sample_name}_initial_alignment.sam 
        samtools view -@ 16 -b ${sample_name}_initial_alignment.sam > ${sample_name}_initial_alignment.bam
        samtools sort -@ 16 ${sample_name}_initial_alignment.bam > ${sample_name}_initial_alignment_sorted.bam
        mv ${sample_name}_initial_alignment_sorted.bam ${sample_name}_initial_alignment.bam
        rm ${sample_name}_initial_alignment.sam
        samtools sort -@ 16 -n -o ${sample_name}_namesort.bam ${sample_name}_initial_alignment.bam
        samtools fixmate -m ${sample_name}_namesort.bam ${sample_name}_fixmate.bam
        samtools sort -@ 16  -o ${sample_name}_positionsort.bam ${sample_name}_fixmate.bam
        samtools markdup  -r ${sample_name}_positionsort.bam ${sample_name}_markdup.bam
        mv ${sample_name}_markdup.bam ${sample_name}_initial_alignment.bam
        samtools index ${sample_name}_initial_alignment.bam
        samtools faidx ${reference}
        rm ${sample_name}_positionsort.bam
        rm ${sample_name}_fixmate.bam
        rm ${sample_name}_namesort.bam
        """
    }

process GET_TARGET_INFORMATION {
    cache 'lenient'
    publishDir "${params.outdir}/${sample_name}/"
    input:
        tuple val(sample_name), path(R1), path(R2), path(metadata_fn)
        val attp_reg
        val attp_prime
    output:
        path "target_info.csv"
    
    script:
    """
    get_target_info.py --metadata ${metadata_fn} --attp_reg ${attp_reg} --attp_prime ${attp_prime}
    """
}

process EXTRACT_TARGET_READS {
    cache 'lenient'
    maxForks 1
    publishDir "${params.outdir}/${sample_name}/"
    input:
        //tuple val(sample_name), path(R1), path(R2), path(metadata_fn)
        path target_info
        val(sample_name)
        tuple path(bam_file), path(bam_file_index)
    output:
        val sample_name, emit: sample_name
        path("${sample_name}_extracted_reads"), emit: extracted_reads_dir
        path("${sample_name}_read_counts_per_site.csv"), emit: read_counts_per_site_file
    
    script:
    """
    extract_target_reads.py --target_info ${target_info} --bam_file ${bam_file} --sample_name ${sample_name}
    """
}

process GENERATE_AMPLICONS {
    cache 'lenient'
    publishDir "${params.outdir}/${sample_name}/"
    input: 
        tuple val(sample_name), path(R1), path(R2), path(metadata)
        path target_info
    output:
        path "amplicons"
    script:
    """
    create_amplicon_files.py --target_info ${target_info}
    """
}

process ALIGN_TARGET_READS {
    cache 'lenient'
    memory 8.GB
    maxForks 1
    publishDir "${params.outdir}/${sample_name}/"
    input:
        // tuple val(sample_name), path(R1), path(R2), path(metadata)
        path target_info
        val(sample_name)
        path(fastq_dir)
        path amplicon_dir
    output:
        val(sample_name), emit: sample_name
        path("${sample_name}_alignments"), emit: probe_read_alignments
    script:
    """
    align_extracted_reads.py --target_info ${target_info} --fastq_dir ${fastq_dir} --amplicon_dir ${amplicon_dir} --sample_name ${sample_name}
    """
}

process MEASURE_INTEGRATION {
    cache 'lenient'
    publishDir "${params.outdir}/${sample_name}", pattern: 'integration_stats.csv'
    input:
        // tuple val(sample_name), path(R1), path(R2), path(metadata)
        path target_info
        val(sample_name)
        path(alignment_dir)
    output:
        val(sample_name), emit: sample_name
        path("${sample_name}_integration_stats.csv"), emit: integration_stats_file
        path("${sample_name}_attL_extracted_reads"), emit: attL_extracted_reads_dir
        path("${sample_name}_attR_extracted_reads"), emit: attR_extracted_reads_dir
        path("${sample_name}_beacon_extracted_reads"), emit: beacon_extracted_reads_dir
    
    script:
    """
    compute_integration_percentage.py --target_info ${target_info} --alignment_dir ${alignment_dir} --sample_name ${sample_name}
    """
}

process GATHER_QC_INFO {
    cache 'lenient'
    publishDir "${params.outdir}/${sample_name}/qc"
    input:
        // tuple val(sample_name), path(R1), path(R2), path(metadata)
        tuple val(sample_name), path(json_file)
        path(fastq_dir)
    output:
        val(sample_name), emit: sample_name
        path "${sample_name}_qc_summary.csv", emit: qc_summary_file
    script:
    """
    gather_qc_stats.py --json_file ${json_file} --fastq_dir ${fastq_dir} --sample_name ${sample_name}
    """
}

process RUN_CS2_FOR_INDELS {
    cache 'lenient'
    maxForks 1
    publishDir "${params.outdir}/${sample_name}/cs2_output"
    input:
        // tuple val(sample_name), path(R1), path(R2), path(metadata)
        val(sample_name)
        path(attL_extracted_reads_dir)
        path(attR_extracted_reads_dir)
        path(beacon_extracted_reads_dir)
        path amplicon_dir
        path target_info
    output:
        val(sample_name), emit: sample_name
        path("${sample_name}_cs2_attL"), emit: cs2_attL_dir
        path("${sample_name}_cs2_attR"), emit: cs2_attR_dir
        path("${sample_name}_cs2_beacon"), emit: cs2_beacon_dir
   
    script:
    """
    run_cs2.py --amplicon_dir ${amplicon_dir} --target_info ${target_info} --sample_name ${sample_name}
    """
}

process GET_INDEL_INFO_FROM_CS2_OUTPUT {
    cache 'lenient'
    maxForks 1
    publishDir "${params.outdir}/${sample_name}/indel_info"
    input:
        //tuple val(sample_name), path(R1), path(R2), path(metadata)
        val(sample_name)
        path(cs2_attL)
        path(cs2_attR)
        path(cs2_beacon)
    output:
        val(sample_name), emit: sample_name
        path("${sample_name}_attL_indel_table.csv"), emit: attL_indel_table
        path("${sample_name}_attR_indel_table.csv"), emit: attR_indel_table
   
    script:
    """
    get_indel_info.py --cs2_directory ${cs2_attL} --output_fn ${sample_name}_attL_indel_table.csv
    get_indel_info.py --cs2_directory ${cs2_attR} --output_fn ${sample_name}_attR_indel_table.csv
    """
}

process COMBINE_INTEGRATION_AND_INDEL_INFO {
    cache 'lenient'
    publishDir "${params.outdir}/${sample_name}"
    input:
        // tuple val(sample_name), path(R1), path(R2), path(metadata)
        // tuple val(sample_name), path(integration_table), path(indel_attR_table), path(indel_attL_table)
        val sample_name
        path integration_table
        path indel_attR_table
        path indel_attL_table
    output:
        val(sample_name), emit: sample_name
        path("${sample_name}_integration_and_indel_stats.csv"), emit: integration_and_indel_stats
    
    script:
    """
    combine_integration_and_indel_stats.py --attL_indel_table ${indel_attL_table} --attR_indel_table ${indel_attR_table} --integration_table ${integration_table} --output_fn ${sample_name}_integration_and_indel_stats.csv 
    """
}

process GENERATE_REPORT {  
    cache 'lenient'
    publishDir "${params.outdir}"
    input:
        path project_config_file
        path integration_stats_files
        path attL_indel_table_files
        path attR_indel_table_files
        path read_counts_per_site_files
        path qc_summary_files
        path extracted_reads_dirs
        val collapse_condition
        val project_name
    output:
        path "*.xlsx", emit: excel_output

    script:
        """
        collate_results.py --project_config_file ${project_config_file} --integration_stats_files ${integration_stats_files} --attL_indel_table_files ${attL_indel_table_files} --attR_indel_table_files ${attR_indel_table_files} --read_counts_per_site_files ${read_counts_per_site_files} --qc_summary_files ${qc_summary_files} --extracted_reads_dirs ${extracted_reads_dirs} --collapse_condition ${collapse_condition} --project_name ${project_name}
        """
    }

process CREATE_PLOTS {  
    cache 'lenient'
    publishDir "${params.outdir}"
    input:
        path excel_file
    output:
        path "*.png"

    script:
        """
        plots.py --excel_report  ${excel_file}
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

// process RUN_MANTA {
//     conda '/home/ubuntu/mambaforge/envs/manta'
//     input:
//         tuple path(bam_file), path(bam_file_index)
//         path reference
//         path reference_index
//     output:
//         tuple path("results/variants/candidateSV.vcf.gz"), path("results/variants/candidateSV.vcf.gz.tbi")
//         // path 'results/evidence/evidence_0.initial_alignment.bam'
//     script:
//     """
//     configManta.py --bam ${bam_file} --referenceFasta ${reference} --runDir . --generateEvidenceBam --exome
//     ./runWorkflow.py
//     """
// }

// process PARSE_MANTA_OUTPUT {
//     input:
//         tuple path(vcf_file), path(vcf_file_index)
//     output:
//         path 'translocation_info.csv'
//     script:
//     """
//     parse_manta_output.py --vcf_file ${vcf_file} --output_csv translocation_info.csv
//     """
// }


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
         project name : ${params.project_name}
         input from   : ${params.samplesheet}
         output to    : ${params.outdir}

         other parameters:
         reference genome: ${params.reference}
         attP, left side: ${params.ATTP_REG}
         attP, right side: ${params.ATTP_PRIME}
         deduplication method: ${params.deduplication_method}
         collapse condition: ${params.collapse_condition}
         --
         run as       : ${workflow.commandLine}
         started at   : ${workflow.start}
         config files : ${workflow.configFiles}
         container    : ${workflow.containerEngine}:${workflow.container}
         """
         .stripIndent()

    Channel.fromPath(params.samplesheet)
    .splitCsv(header: true, sep: ',')
    .map { row -> 
        def r1 = "${launchDir}/${row.fastq_dir}/*R1*.fastq.gz"
        def r2 = "${launchDir}/${row.fastq_dir}/*R2*.fastq.gz"
        tuple(
            row.sample_name, 
            file(r1), 
            file(r2), 
            file(row.probe_list)
        )
    }
    .set { input_ch }

    def reference_absolute_path = "${launchDir}/${params.reference}"

    trimmed_and_merged_fastq = ADAPTER_AND_POLY_G_TRIM(input_ch)

    minimap2_reference_index = CREATE_MINIMAP2_REFERENCE_INDEX(reference_absolute_path)

    initial_alignment = ALIGN_READS(reference_absolute_path, minimap2_reference_index, trimmed_and_merged_fastq.trimmed_fastq, params.deduplication_method)

    probe_information = GET_TARGET_INFORMATION(input_ch, params.ATTP_REG, params.ATTP_PRIME)

    fastq_dir = EXTRACT_TARGET_READS(probe_information, initial_alignment.sample_name, initial_alignment.initial_alignment_bam)

    amplicons_dir = GENERATE_AMPLICONS(input_ch,probe_information)

    alignment_dir = ALIGN_TARGET_READS(probe_information, fastq_dir.sample_name, fastq_dir.extracted_reads_dir, amplicons_dir)

    att_sequence_dirs = MEASURE_INTEGRATION(probe_information,alignment_dir.sample_name, alignment_dir.probe_read_alignments)

    qc_summary = GATHER_QC_INFO(trimmed_and_merged_fastq.fastp_stats, fastq_dir.extracted_reads_dir)

    cs2_attL_attR_dirs = RUN_CS2_FOR_INDELS(att_sequence_dirs.sample_name, att_sequence_dirs.attL_extracted_reads_dir, att_sequence_dirs.attR_extracted_reads_dir, att_sequence_dirs.beacon_extracted_reads_dir, amplicons_dir, probe_information)

    indel_tables = GET_INDEL_INFO_FROM_CS2_OUTPUT(cs2_attL_attR_dirs.sample_name, cs2_attL_attR_dirs.cs2_attL_dir, cs2_attL_attR_dirs.cs2_attR_dir, cs2_attL_attR_dirs.cs2_beacon_dir)

    COMBINE_INTEGRATION_AND_INDEL_INFO(att_sequence_dirs.sample_name, att_sequence_dirs.integration_stats_file, indel_tables.attL_indel_table, indel_tables.attR_indel_table)

    def samplesheet_absolute_path = "${launchDir}/${params.samplesheet}"
    
    report_excel_file = GENERATE_REPORT(samplesheet_absolute_path,att_sequence_dirs.integration_stats_file.collect(), indel_tables.attL_indel_table.collect(), indel_tables.attR_indel_table.collect(), fastq_dir.read_counts_per_site_file.collect(), qc_summary.qc_summary_file.collect(), fastq_dir.extracted_reads_dir.collect(), params.collapse_condition, params.project_name)
    
    MULTIQC(trimmed_and_merged_fastq.fastp_stats
        .flatten()  // Flatten the list
        .filter { it.toString().endsWith('.json') }  // Filter out only the paths ending with .json
        .collect()
        .ifEmpty([]))

    CREATE_PLOTS(report_excel_file.excel_output)
    
}

workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()
}
