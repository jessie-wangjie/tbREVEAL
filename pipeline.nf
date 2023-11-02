nextflow.enable.dsl=2

params.reference = "/data/references/hg38.fa"
params.cargo = '/data/references/AAVG107.fa'
params.outdir = ""
params.ATTP_REG = 'GTGGTTTGTCTGGTCAACCACCGCGGT'
params.ATTP_PRIME = 'CTCAGTGGTGTACGGTACAAACCCA'
params.deduplication_method = 'STRICT'
params.collapse_condition = 'Complete'
params.initial_mapper = 'minimap2'
params.project_name = ''
params.samplesheet = ''
params.reference_index_location = ''
params.lmpcr_mode = false

process ADAPTER_AND_POLY_G_TRIM {
    cache 'lenient'
    maxForks 1
    // publishDir "${params.outdir}/${sample_name}/qc/", pattern: '*fastp.json'
    publishDir "${params.outdir}/raw_fastq/${sample_name}/", pattern: '*.fastq.gz'
    publishDir "${params.outdir}/trimmed_fastq/${sample_name}/", pattern: '*.fastq.bgz'
    publishDir "${params.outdir}/input_probe_sheet/${sample_name}/", pattern: '*.csv'

    input:
        tuple val(sample_name), path(R1), path(R2), path(metadata), val(group)
        val lmpcr_mode
    output:
        path "${R1}"
        path "${R2}"
        path "${metadata}"
        tuple val(sample_name), path("${sample_name}_trimmed_R1.fastq.gz"), path("${sample_name}_trimmed_R2.fastq.gz"), emit: trimmed_nonmerged_fastq
        tuple val(sample_name), val(group), path("${sample_name}_trimmed.fastq.bgz"), emit: trimmed_fastq
        tuple val(sample_name), path("${sample_name}_fastp.json"), emit: fastp_stats

    script:
        if (lmpcr_mode == false) {
        """
        fastp -m -c --include_unmerged --low_complexity_filter --cut_right --overlap_len_require 10 -i ${R1} -I ${R2} --merged_out ${sample_name}_trimmed.fastq.gz -w 16 -g -j ${sample_name}_fastp.json
        gunzip -c ${sample_name}_trimmed.fastq.gz | bgzip -@ 8 > ${sample_name}_trimmed.fastq.bgz

        fastp --low_complexity_filter --cut_right -i ${R1} -I ${R2} -o ${sample_name}_trimmed_R1.fastq.gz -O ${sample_name}_trimmed_R2.fastq.gz -w 16 -g
        """
        } else {
        """
        fastp -m -c --include_unmerged --low_complexity_filter --cut_right --overlap_len_require 10 -i ${R1} -I ${R2} --merged_out ${sample_name}_trimmed.fastq.gz -w 16 -g -j ${sample_name}_fastp.json 

        fastp --low_complexity_filter --cut_right -i ${R1} -I ${R2} -o ${sample_name}_trimmed_R1.fastq.gz -O ${sample_name}_trimmed_R2.fastq.gz -w 16 -g

        gunzip ${sample_name}_trimmed.fastq.gz

        AmpUMI Process --fastq ${sample_name}_trimmed.fastq --fastq_out ${sample_name}_trimmed2.fastq --umi_regex "^IIIIIIIIIII"

        mv ${sample_name}_trimmed2.fastq ${sample_name}_trimmed.fastq

        gzip ${sample_name}_trimmed.fastq

        gunzip -c ${sample_name}_trimmed.fastq.gz | bgzip -@ 8 > ${sample_name}_trimmed.fastq.bgz
        """
        }
        
}

process CREATE_REFERENCE_INDEX {
    cache 'lenient'

    input:
        path reference
        val mapper

    output:
        path ['*.fa.*', '*.mmi'], emit: reference_index_files

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
    maxForks 1
    publishDir "${params.outdir}/initial_alignments/${sample_name}/", pattern: '*_initial_alignment.bam*'
    publishDir "${params.outdir}/deduped_alignments/${sample_name}/", pattern: '*_deduped_alignment.bam*'

    input:
        path reference
        path reference_index
        tuple val(sample_name), val(group), path(fastq)
        val dedup_method
        val initial_mapper

    output:
        val(sample_name), emit: sample_name
        tuple val(group), path("${sample_name}_initial_alignment.bam"), path("${sample_name}_initial_alignment.bam.bai")
        tuple val(group), path("${sample_name}_deduped_alignment.bam"), path("${sample_name}_deduped_alignment.bam.bai"), emit: deduped_alignment_bam
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

    if (dedup_method == "STRICT") {
        """
        $alignment_command
        samtools view -@ 16 -b ${sample_name}_initial_alignment.sam > ${sample_name}_initial_alignment.bam
        samtools sort -@ 16 ${sample_name}_initial_alignment.bam > ${sample_name}_initial_alignment_sorted.bam
        mv ${sample_name}_initial_alignment_sorted.bam ${sample_name}_initial_alignment.bam
        rm ${sample_name}_initial_alignment.sam
        samtools index ${sample_name}_initial_alignment.bam
        samtools sort -@ 16 -n -o ${sample_name}_namesort.bam ${sample_name}_initial_alignment.bam
        samtools fixmate -m ${sample_name}_namesort.bam ${sample_name}_fixmate.bam
        samtools sort -@ 16  -o ${sample_name}_positionsort.bam ${sample_name}_fixmate.bam
        samtools markdup --barcode-name -r ${sample_name}_positionsort.bam ${sample_name}_markdup.bam  
        mv ${sample_name}_markdup.bam ${sample_name}_deduped_alignment.bam
        samtools index ${sample_name}_deduped_alignment.bam
        samtools faidx ${reference}
        rm ${sample_name}_positionsort.bam
        rm ${sample_name}_fixmate.bam
        rm ${sample_name}_namesort.bam
        """
    } else if (dedup_method == "LOOSE") {
        """
        $alignment_command
        samtools view -@ 16 -b ${sample_name}_initial_alignment.sam > ${sample_name}_initial_alignment.bam
        samtools sort -@ 16 ${sample_name}_initial_alignment.bam > ${sample_name}_initial_alignment_sorted.bam
        mv ${sample_name}_initial_alignment_sorted.bam ${sample_name}_initial_alignment.bam
        rm ${sample_name}_initial_alignment.sam
        samtools index ${sample_name}_initial_alignment.bam
        samtools sort -@ 16 -n -o ${sample_name}_namesort.bam ${sample_name}_initial_alignment.bam
        samtools fixmate -m ${sample_name}_namesort.bam ${sample_name}_fixmate.bam
        samtools sort -@ 16  -o ${sample_name}_positionsort.bam ${sample_name}_fixmate.bam
        samtools markdup  -r ${sample_name}_positionsort.bam ${sample_name}_markdup.bam
        mv ${sample_name}_markdup.bam ${sample_name}_deduped_alignment.bam
        samtools index ${sample_name}_deduped_alignment.bam
        samtools faidx ${reference}
        rm ${sample_name}_positionsort.bam
        rm ${sample_name}_fixmate.bam
        rm ${sample_name}_namesort.bam
        """
    } else if (dedup_method == "LMPCR") {
        """
        $alignment_command
        samtools view -@ 16 -b ${sample_name}_initial_alignment.sam > ${sample_name}_initial_alignment.bam
        samtools sort -@ 16 ${sample_name}_initial_alignment.bam > ${sample_name}_initial_alignment_sorted.bam
        mv ${sample_name}_initial_alignment_sorted.bam ${sample_name}_initial_alignment.bam
        rm ${sample_name}_initial_alignment.sam
        samtools index ${sample_name}_initial_alignment.bam
        cp ${sample_name}_initial_alignment.bam ${sample_name}_deduped_alignment.bam
        samtools index ${sample_name}_deduped_alignment.bam
        samtools faidx ${reference}
        """
    } 
}

process GET_TARGET_INFORMATION {
    cache 'lenient'
    publishDir "${params.outdir}/full_probe_info/${sample_name}/"
    input:
        tuple val(sample_name), path(R1), path(R2), path(metadata_fn), val(group)
        path reference
        path cargo_ref
        val attp_reg
        val attp_prime
    output:
        path "target_info.csv"
    
    script:
    """
    get_target_info.py --metadata ${metadata_fn} --attp_reg ${attp_reg} --attp_prime ${attp_prime} --reference ${reference} --cargo ${cargo_ref}
    """
}

process EXTRACT_TARGET_READS {
    cache 'lenient'
    maxForks 1
    publishDir "${params.outdir}/probe_specific_reads/${sample_name}/"
    input:
        path target_info
        val(sample_name)
        tuple val(group), path(bam_file), path(bam_file_index)
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
    publishDir "${params.outdir}/generated_amplicons/${sample_name}/"
    input: 
        path target_info
        val(sample_name)
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
    publishDir "${params.outdir}/read_to_probe_alignment/"
    input:
        path target_info
        val sample_name 
        path fastq_dir 
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
    publishDir "${params.outdir}/integration_stats_tables/${sample_name}", pattern: 'integration_stats.csv'
    input:
        path target_info
        val sample_name 
        path alignment_dir 
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
    publishDir "${params.outdir}/qc_info/${sample_name}/"
    input:
        tuple val(sample_name), path(json_file)
        path fastq_dir 
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
    // publishDir "${params.outdir}/${sample_name}/cs2_output",pattern:"*cs2_*"
    publishDir "${params.outdir}/alignment_visualizations/", pattern:"*alignments*"

    input:
        val(sample_name)
        path attL_extracted_reads_dir 
        path attR_extracted_reads_dir 
        path beacon_extracted_reads_dir 
        path amplicon_dir
        path target_info
    output:
        val(sample_name), emit: sample_name
        path("${sample_name}_cs2_attL"), emit: cs2_attL_dir
        path("${sample_name}_cs2_attR"), emit: cs2_attR_dir
        path("${sample_name}_cs2_beacon"), emit: cs2_beacon_dir

        path("${sample_name}_attL_alignments"), optional: true, emit: cs2_attL_alignments
        path("${sample_name}_attR_alignments"), optional: true, emit: cs2_attR_alignments
        path("${sample_name}_beacon_alignments"), optional: true, emit: cs2_beacon_alignments
   
    script:
    """
    run_cs2.py --amplicon_dir ${amplicon_dir} --target_info ${target_info} --sample_name ${sample_name}
    """
}

process GET_INDEL_INFO_FROM_CS2_OUTPUT {
    cache 'lenient'
    maxForks 1
    publishDir "${params.outdir}/indel_info/${sample_name}/"
    input:
        val sample_name 
        path cs2_attL 
        path cs2_attR 
        path cs2_beacon 
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
    // publishDir "${params.outdir}/${sample_name}"
    input:
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

process MERGE_ALIGNMENTS {
    cache 'lenient'
    publishDir "${params.outdir}/deduped_alignments_merged_by_group/"
    input:
        tuple val(group), path(bam_files), path(bai_files)

    output:
        path("${group}_samples_alignment.bam*")

    script:
        """
        samtools merge ${group}_samples_alignment.bam ${bam_files.join(' ')}
        samtools index ${group}_samples_alignment.bam
        """
}

process RUN_SV_CALLING {
    maxForks 1
    cache 'lenient'
    publishDir "${params.outdir}/sv_output/raw_vcf/"
    input:
        tuple val(sample_name), path(R1), path(R2)
        path reference
        path reference_index
    output:
        tuple val(sample_name), path("*dysgu_svs.vcf"), emit: sv_output
    script:
    """
    minimap2 -ax sr ${reference} ${R1} ${R2} | samtools view -b - | samtools sort -@ 16 - > ${sample_name}_unmerged_reads_aligned_to_reference.bam
    samtools index ${sample_name}_unmerged_reads_aligned_to_reference.bam
    dysgu run --max-cov 99999 --overwrite ${reference} temp_dir ${sample_name}_unmerged_reads_aligned_to_reference.bam > ${sample_name}_dysgu_svs.vcf
    """
}

process PARSE_SV_CALLING_OUTPUT {
    cache 'lenient'
    publishDir "${params.outdir}/sv_output/translocation_analysis/"
    input:
        tuple val(sample_name), path(vcf_file)
    output:
        path '*translocation_info.csv'
    script:
    """
    parse_manta_output.py --vcf ${vcf_file} --output_csv ${sample_name}_translocation_info.csv
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

         other parameters:
         reference genome: ${params.reference}
         initial mapper: ${params.initial_mapper}
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
            file(row.probe_list),
            row.group
        )
    }
    .set { input_ch }

    Channel.fromPath("${params.reference_index_location}/*.{fa,fasta}.*")
    .collect()
    .set { reference_index_ch }

    def reference_absolute_path = "${launchDir}/${params.reference}"
    def cargo_absolute_path = "${launchDir}/${params.cargo}"

    if (params.reference_index_location == '') {
        create_reference_index = CREATE_REFERENCE_INDEX(reference_absolute_path, params.initial_mapper)
        reference_index_ch = create_reference_index.reference_index_files
    }

    trimmed_and_merged_fastq = ADAPTER_AND_POLY_G_TRIM(input_ch, params.lmpcr_mode)

    initial_alignment = ALIGN_READS(reference_absolute_path, reference_index_ch, trimmed_and_merged_fastq.trimmed_fastq, params.deduplication_method,params.initial_mapper)

    probe_information = GET_TARGET_INFORMATION(input_ch, reference_absolute_path, cargo_absolute_path, params.ATTP_REG, params.ATTP_PRIME)

    fastq_dir = EXTRACT_TARGET_READS(probe_information, initial_alignment.sample_name, initial_alignment.deduped_alignment_bam)

    amplicons_dir = GENERATE_AMPLICONS(probe_information,initial_alignment.sample_name)

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


    initial_alignment.deduped_alignment_bam
        .groupTuple()
        .set { grouped_bam_files }
    

    MERGE_ALIGNMENTS(grouped_bam_files)

    sv_calling = RUN_SV_CALLING(trimmed_and_merged_fastq.trimmed_nonmerged_fastq, reference_absolute_path, initial_alignment.reference_fasta_fai)

    PARSE_SV_CALLING_OUTPUT(sv_calling.sv_output)
    
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
