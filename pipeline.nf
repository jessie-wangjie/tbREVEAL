nextflow.enable.dsl=2

params.reference = "/data/references/hg38.fa"
params.outdir = ""
params.ATTP_REG = 'GTGGTTTGTCTGGTCAACCACCGCGGT'
params.ATTP_PRIME = 'CTCAGTGGTGTACGGTACAAACCCA'
params.collapse_condition = 'Complete'
params.initial_mapper = 'minimap2'
params.project_name = ''
params.samplesheet = ''
params.reference_index_location = ''
params.lmpcr_mode = false
params.useBam = false
params.umi_in_header = false
params.umi_loc = 'per_read'
params.umi_length = 5
params.other_fastp_params = ''
params.notebook_template = "${workflow.projectDir}/bin/report_generation.ipynb"
params.bam2html_path = "${workflow.projectDir}/bin/utils/bam2html.py"

process ADAPTER_AND_POLY_G_TRIM {
    cache 'lenient'
    publishDir "${params.outdir}/raw_fastq/${sample_name}/", pattern: '*.fastq.gz'
    publishDir "${params.outdir}/trimmed_fastq/${sample_name}/", pattern: '*.fastq.bgz'
    publishDir "${params.outdir}/input_probe_sheet/${sample_name}/", pattern: '*.csv'

    input:
        tuple val(sample_name), path(R1), path(R2), path(metadata), path(cargo_ref),val(group)
        val umi_in_header
        val umi_loc
        val umi_length
        val other_fastp_params
    output:
        path "${R1}"
        path "${R2}"
        path "${metadata}"
        tuple val(sample_name), val(group), path("${sample_name}_trimmed.fastq.bgz"), emit: trimmed_fastq
        tuple val(sample_name), val(group), path("${sample_name}_fastp.json"), emit: fastp_stats

    script:
        if (umi_in_header == false) {
            """
            fastp -m -c --dont_eval_duplication --disable_adapter_trimming --low_complexity_filter ${other_fastp_params} --overlap_len_require 10 -i ${R1} -I ${R2} --merged_out ${sample_name}_trimmed.fastq.gz -w 16 -g -j ${sample_name}_fastp.json -U --umi_loc=${umi_loc} --umi_len=${umi_length} --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

            gunzip -c ${sample_name}_trimmed.fastq.gz | bgzip -@ 8 > ${sample_name}_trimmed.fastq.bgz

            # fastp --low_complexity_filter --dont_eval_duplication ${other_fastp_params} -i ${R1} -I ${R2} -o ${sample_name}_trimmed_R1.fastq.gz -O ${sample_name}_trimmed_R2.fastq.gz -w 16 -U --umi_loc=${umi_loc} --umi_len=${umi_length}
            """
        } else {
            """
            fastp -m -c --dont_eval_duplication --disable_adapter_trimming --low_complexity_filter ${other_fastp_params} --overlap_len_require 10 -i ${R1} -I ${R2} --merged_out ${sample_name}_trimmed.fastq.gz -w 16 -g -j ${sample_name}_fastp.json --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

            gunzip -c ${sample_name}_trimmed.fastq.gz | bgzip -@ 8 > ${sample_name}_trimmed.fastq.bgz

            # fastp --low_complexity_filter --dont_eval_duplication ${other_fastp_params} -i ${R1} -I ${R2} -o ${sample_name}_trimmed_R1.fastq.gz -O ${sample_name}_trimmed_R2.fastq.gz -w 16
            """
        }
    }

process CREATE_REFERENCE_INDEX {
    cache 'lenient'

    input:
        path reference
        val mapper

    output:
        path ['*.fa.*', '*.fasta.*', '*.mmi'], emit: reference_index_files

    script:
    if (mapper == "minimap2") {
        """
        minimap2 -x sr -d reference_index.mmi ${reference}
        """
    } else if (mapper == "bwa") {
        """
        bwa index ${reference}
        """
    } else {
        error "Unsupported mapper: $mapper"
    }
}

process ALIGN_READS {
    cache 'lenient'
    publishDir "${params.outdir}/initial_alignments/${sample_name}/", pattern: '*_initial_alignment.bam*'
    publishDir "${params.outdir}/deduped_alignments/${sample_name}/", pattern: '*_deduped_alignment.bam*'

    input:
        path reference
        path reference_index
        tuple val(sample_name), val(group), path(fastq)
        val initial_mapper
        val umi_deduplication

    output:
        tuple val(sample_name), val(group), path("${sample_name}_initial_alignment.bam"), path("${sample_name}_initial_alignment.bam.bai"), emit: original_alignment_bam
        tuple val(sample_name), val(group), path("${sample_name}_deduped_alignment.bam"), path("${sample_name}_deduped_alignment.bam.bai"), emit: deduped_alignment_bam
        path "${reference}.fai", emit: reference_fasta_fai

    script:
    def alignment_command = ""
    if (initial_mapper == "minimap2") {
        alignment_command = "minimap2 -ax sr -t 16 ${reference_index} ${fastq} > ${sample_name}_initial_alignment.sam"
    } else if (initial_mapper == "bwa") {
        alignment_command = "bwa mem -t 16 ${reference} ${fastq} > ${sample_name}_initial_alignment.sam"
    } else {
        error "Unsupported initial_mapper: $initial_mapper"
    }
    if (umi_deduplication == true) {
        """
        $alignment_command
        samtools view -@ 16 -b ${sample_name}_initial_alignment.sam > ${sample_name}_initial_alignment.bam
        samtools sort -@ 16 ${sample_name}_initial_alignment.bam > ${sample_name}_initial_alignment_sorted.bam
        mv ${sample_name}_initial_alignment_sorted.bam ${sample_name}_initial_alignment.bam
        samtools index ${sample_name}_initial_alignment.bam
        umi_tools dedup -I ${sample_name}_initial_alignment.bam --paired --umi-separator ":" -S ${sample_name}_deduped_alignment.bam --method unique
        rm ${sample_name}_initial_alignment.sam
        samtools index ${sample_name}_deduped_alignment.bam
        samtools faidx ${reference}
        """
    } else {
        """
        $alignment_command
        samtools view -@ 16 -b ${sample_name}_initial_alignment.sam > ${sample_name}_initial_alignment.bam
        samtools sort -@ 16 ${sample_name}_initial_alignment.bam > ${sample_name}_initial_alignment_sorted.bam
        mv ${sample_name}_initial_alignment_sorted.bam ${sample_name}_initial_alignment.bam
        samtools index ${sample_name}_initial_alignment.bam
        cp ${sample_name}_initial_alignment.bam ${sample_name}_deduped_alignment.bam
        rm ${sample_name}_initial_alignment.sam
        samtools index ${sample_name}_deduped_alignment.bam
        samtools faidx ${reference}
        """
    }
}

process GET_TARGET_INFORMATION {
    cache 'lenient'
    publishDir "${params.outdir}/full_probe_info/${sample_name}/"
    input:
        tuple val(sample_name), path(R1), path(R2), path(metadata_fn), path(cargo_ref),val(group)
        path reference
        val attp_reg
        val attp_prime
    output:
        tuple val(sample_name), val(group), path("${sample_name}_target_info.csv")

    script:
    """
    get_target_info.py --metadata ${metadata_fn} --attp_reg ${attp_reg} --attp_prime ${attp_prime} --reference ${reference} --cargo ${cargo_ref} --sample_name ${sample_name}
    """
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
        path bam2html_path
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
        collate_results.py --project_config_file ${project_config_file} --integration_stats_files ${integration_stats_files} --read_counts_per_site_files ${read_counts_per_site_files} --qc_summary_files ${qc_summary_files} --extracted_reads_dirs ${extracted_reads} --collapse_condition ${collapse_condition} --project_name ${project_name}
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

process CREATE_PYTHON_NOTEBOOK_REPORT {
    cache 'lenient'
    publishDir "${params.outdir}"

    input:
        path excel_file
        path notebook_template
    output:
        path 'report.html'
    script:
    """
    papermill ${notebook_template} report.ipynb -p results_file ${excel_file}
    jupyter nbconvert --to html --no-input report.ipynb
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
         project name : ${params.project_name}
         input from   : ${params.samplesheet}
         output to    : ${params.outdir}
         workflow directory : ${workflow.projectDir}

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

    Channel.fromPath("${params.reference_index_location}/*.{fa,fasta,fna,mmi}.*")
    .collect()
    .set { reference_index_ch }


    def reference_absolute_path = "${launchDir}/${params.reference}"

    if (params.reference_index_location == '') {
        create_reference_index = CREATE_REFERENCE_INDEX(reference_absolute_path, params.initial_mapper)
        reference_index_ch = create_reference_index.reference_index_files
    }

    if (params.useBam) {
    // Logic to handle BAM file input
        Channel.fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def bam = "${launchDir}/${row.fastq_dir}/*.bam"
            def bam_index = "${launchDir}/${row.fastq_dir}/*.bam.bai"
            tuple(
                row.sample_name,
                file(bam),
                file(bam_index),
                file(row.probe_list),
                row.group
        )
    }
    .set { bam_input_ch }
    } else {
        // Existing logic to handle FASTQ file input
        Channel.fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def r1 = "${launchDir}/${row.fastq_dir}/*R1*.fastq.gz"
            def r2 = "${launchDir}/${row.fastq_dir}/*R2*.fastq.gz"
            tuple(
                row.sample_name,
                file(r1),
                file(r2),
                file(row.probe_list),
                file(row.cargo),
                row.group
            )
        }
        .set { input_ch }
    }



    if (params.useBam){
        sample_name = bam_input_ch.map {
           it[0]
        }
        bam_file = bam_input_ch.map {
           it[1]
        }
        bam_index = bam_input_ch.map {
           it[2]
        }
        group = bam_input_ch.map {
           it[4]
        }

        bam_tuple_ch = bam_input_ch.map { items ->
            return tuple(items[0], items[4], items[1], items[2])
        }

        probe_information = GET_TARGET_INFORMATION(bam_input_ch, reference_absolute_path, cargo_absolute_path, params.ATTP_REG, params.ATTP_PRIME)
        fastq_dir = EXTRACT_TARGET_READS(probe_information, sample_name, bam_tuple_ch)
        amplicons_dir = GENERATE_AMPLICONS(probe_information, sample_name)
        alignment_dir = ALIGN_TARGET_READS(probe_information, fastq_dir.sample_name, fastq_dir.extracted_reads_dir, amplicons_dir)
        att_sequence_dirs = MEASURE_INTEGRATION(probe_information,alignment_dir.sample_name, alignment_dir.probe_read_alignments)
        cs2_attL_attR_dirs = ALIGNMENT_VISUALIZATION(att_sequence_dirs.sample_name, alignment_dir.probe_read_alignments, att_sequence_dirs.attL_extracted_reads_dir, att_sequence_dirs.attR_extracted_reads_dir, att_sequence_dirs.beacon_extracted_reads_dir, att_sequence_dirs.wt_extracted_reads_dir, amplicons_dir, probe_information)

        // def samplesheet_absolute_path = "${launchDir}/${params.samplesheet}"
        // report_excel_file = GENERATE_REPORT(samplesheet_absolute_path,att_sequence_dirs.integration_stats_file.collect(), indel_tables.attL_indel_table.collect(), indel_tables.attR_indel_table.collect(), fastq_dir.read_counts_per_site_file.collect(), qc_summary.qc_summary_file.collect(), fastq_dir.extracted_reads_dir.collect(), params.collapse_condition, params.project_name)
        // CREATE_PLOTS(report_excel_file.excel_output)

    } else {
        // run when fastq input

        // *** GET PROBE INFO ***
        probe_information = GET_TARGET_INFORMATION(input_ch, reference_absolute_path, params.ATTP_REG, params.ATTP_PRIME)

        // *** CLEAN READS ***
        trimmed_and_merged_fastq = ADAPTER_AND_POLY_G_TRIM(input_ch, params.umi_in_header, params.umi_loc, params.umi_length,params.other_fastp_params)

        // *** ALIGN READS ***
        initial_alignment = ALIGN_READS(reference_absolute_path, reference_index_ch, trimmed_and_merged_fastq.trimmed_fastq, params.initial_mapper, params.umi_deduplication)

        // ** EXTRACT PROBE READS **
        // by: [0,1] combines the channels using sample_name,group as key
        probe_information
            .combine(initial_alignment.deduped_alignment_bam, by: [0,1])
            .set{target_info_and_deduped_alignment_ch}
        extract_target_reads_out = EXTRACT_TARGET_READS(target_info_and_deduped_alignment_ch)

        // ** GENERATE AMPLICONS **
        amplicon_files = GENERATE_AMPLICONS(probe_information)

        // ** ALIGN PROBE READS TO AMPLICONS **
        probe_information
            .combine(extract_target_reads_out.extracted_reads, by: [0,1])
            .combine(amplicon_files, by:[0,1])
            .set{align_target_reads_input_ch}

        align_target_reads_out = ALIGN_TARGET_READS(align_target_reads_input_ch)

        // ** MEASURE INTEGRATION **
        probe_information
            .combine(align_target_reads_out,by:[0,1])
            .set{measure_integration_input_ch}
        measure_integration_out = MEASURE_INTEGRATION(measure_integration_input_ch)

        // ** SUMMARIZE QC STATS **

        trimmed_and_merged_fastq.fastp_stats
            .combine(initial_alignment.original_alignment_bam,by:[0,1])
            .combine(initial_alignment.deduped_alignment_bam,by:[0,1])
            .set{gather_qc_info_input_ch}
        qc_summary = GATHER_QC_INFO(gather_qc_info_input_ch)

        // ** CREATE ALIGNMENT VISUALIZATION **

        probe_information
            .combine(amplicon_files,by:[0,1])
            .combine(align_target_reads_out,by:[0,1])
            .combine(measure_integration_out.edited_reads,by:[0,1])
            .set{alignment_viz_input_ch}
        ALIGNMENT_VISUALIZATION(alignment_viz_input_ch, params.bam2html_path)

        // ** CREATE COLLATED REPORT **

        def samplesheet_absolute_path = "${launchDir}/${params.samplesheet}"

        measure_integration_out.integration_stats_file
            .collect(flat:false)
            .flatMap{ it }
            .map{ tuple -> tuple[2]}
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

        report_excel_file = GENERATE_REPORT(samplesheet_absolute_path,integration_stats_files_ch,read_counts_per_site_files_ch,qc_summary_files_ch,extracted_reads_files_ch, params.collapse_condition, params.project_name)

        // ** MULTIQC REPORT **
        trimmed_and_merged_fastq.fastp_stats
            .flatten()  // Flatten the list
            .filter { it.toString().endsWith('.json') }  // Filter out only the paths ending with .json
            .collect()
            .ifEmpty([])
            .set{multiqc_input_ch}
        MULTIQC(multiqc_input_ch)

        // ** CREATE HTML REPORT **

        CREATE_PYTHON_NOTEBOOK_REPORT(report_excel_file, params.notebook_template)
    }

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
