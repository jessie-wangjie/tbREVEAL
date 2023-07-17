nextflow.enable.dsl=2

params.r1 = "/data/reads/CM_Custom_TB000175c_20230613/4L-550-rep1_L1_ds.228ab3596da04f878fcc010e9dc007b5/4L-550-rep1_S1_L001_R1_001.fastq.gz"
params.r2 = "/data/reads/CM_Custom_TB000175c_20230613/4L-550-rep1_L1_ds.228ab3596da04f878fcc010e9dc007b5/4L-550-rep1_S1_L001_R2_001.fastq.gz"
params.metadata = "/data/CM_Custom_TB000175c_20230613_analysis/20230329-final_panel_up.csv"
params.reference = "/data/references/hg38.fa"
params.outdir = "/data/4L-550-rep1/"
params.sample_name = "test"
params.ATTP_REG = 'GTGGTTTGTCTGGTCAACCACCGCGGT'
params.ATTP_PRIME = 'CTCAGTGGTGTACGGTACAAACCCA'

process ADAPTER_AND_POLY_G_TRIM {
     publishDir "${params.outdir}/${params.sample_name}/qc/", pattern: 'fastp.json'

    input:
        path polyg_trim_input_r1
        path polyg_trim_input_r2

    output:
        tuple path("trimmed_R1.fastq.bgz"), path("trimmed_R2.fastq.bgz")
        path("fastp.json")

    script:
        """
        fastp -i ${polyg_trim_input_r1} -I ${polyg_trim_input_r2} -o trimmed_R1.fastq.gz -O trimmed_R2.fastq.gz -w 32 -g
        gunzip -c trimmed_R1.fastq.gz | bgzip -@ 8 > trimmed_R1.fastq.bgz
        gunzip -c trimmed_R2.fastq.gz | bgzip -@ 8 > trimmed_R2.fastq.bgz
        """
}


process ALIGN_READS {
    publishDir "${params.outdir}/${params.sample_name}/initial_alignment/"

    input:
        path reference
        tuple path(fastq1_file), path(fastq2_file)

    output:
        tuple path("output_sorted.bam"), path("output_sorted.bam.bai")

    script:
        """
        minimap2 -ax sr -t 32 ${reference} ${fastq1_file} ${fastq2_file} \
        | samtools view -@ 8 -b - \
        | samtools sort -@ 8 -o output_sorted.bam -
        samtools index output_sorted.bam
        """
}

process GET_TARGET_INFORMATION {
    publishDir "${params.outdir}/${params.sample_name}/target_info/"
    input:
        path metadata_fn
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
    publishDir "${params.outdir}/${params.sample_name}/"
    input:
        path target_info
        tuple path(bam_file), path(bam_file_index)
        tuple path(fastq1_file), path(fastq2_file)
    output:
        path "extracted_reads"
    
    script:
    """
    extract_target_reads.py --target_info ${target_info} --bam_file ${bam_file}
    """
}

process GENERATE_AMPLICONS {
    publishDir "${params.outdir}/${params.sample_name}/"
    input: 
        path target_info
    output:
        path "amplicons"
    script:
    """
    create_amplicon_files.py --target_info ${target_info}
    """
}

process ALIGN_TARGET_READS {
    publishDir "${params.outdir}/${params.sample_name}/"
    input:
        path target_info
        path fastq_dir
        path amplicon_dir
    output:
        path "alignments"
    script:
    """
    align_extracted_reads.py --target_info ${target_info} --fastq_dir ${fastq_dir} --amplicon_dir ${amplicon_dir}
    """
}

process MEASURE_INTEGRATION {
    publishDir "${params.outdir}/${params.sample_name}/"
    input:
        path target_info
        path alignment_dir
    output:
        path "integration_stats.csv" 
        path "attL_extracted_reads" 
        path "attR_extracted_reads" 
    
    script:
    """
    compute_integration_percentage.py --target_info ${target_info} --alignment_dir ${alignment_dir}
    """
}

process GATHER_QC_INFO {
    publishDir "${params.outdir}/${params.sample_name}/qc"
    input:
        path json_file
        path fastq_dir
    output:
        path 'qc_summary.csv'
        path 'probe_read_counts.csv'
    script:
    """
    gather_qc_stats.py --json_file ${json_file} --fastq_dir ${fastq_dir}
    """
}

process RUN_CS2_FOR_INDELS {
    publishDir "${params.outdir}/${params.sample_name}/cs2_output"
    input:
        path attL_sequences
        path attR_sequences
        path amplicon_dir
        path target_info
    output:
        path 'cs2_attL'
        path 'cs2_attR'
   
    script:
    """
    run_cs2.py --amplicon_dir ${amplicon_dir} --target_info ${target_info}
    """
}

process GET_INDEL_INFO_FROM_CS2_OUTPUT {
    publishDir "${params.outdir}/${params.sample_name}/indel_info"
    input:
        path cs2_attL
        path cs2_attR
    output:
        path "attL_indel_table.csv"
        path "attR_indel_table.csv"
   
    script:
    """
    get_indel_info.py --cs2_directory ${cs2_attL} --output_fn attL_indel_table.csv
    get_indel_info.py --cs2_directory ${cs2_attR} --output_fn attR_indel_table.csv
    """
}

process COMBINE_INTEGRATION_AND_INDEL_INFO {
    publishDir "${params.outdir}/${params.sample_name}"
    input:
        path integration_table
        path indel_attR_table
        path indel_attL_table
    output:
        path "integration_and_indel_stats.csv"
    
    script:
    """
    combine_integration_and_indel_stats.py --attL_indel_table ${indel_attL_table} --attR_indel_table ${indel_attR_table} --integration_table ${integration_table} --output_fn integration_and_indel_stats.csv
    """
}

workflow {
    trim_outputs = ADAPTER_AND_POLY_G_TRIM(params.r1, params.r2)
    align_reads_output = ALIGN_READS(params.reference, trim_outputs[0])
    target_information = GET_TARGET_INFORMATION(params.metadata, params.ATTP_REG, params.ATTP_PRIME)
    fastq_dir = EXTRACT_TARGET_READS(target_information, align_reads_output, trim_outputs[0])
    amplicons_dir = GENERATE_AMPLICONS(target_information)
    alignment_dir = ALIGN_TARGET_READS(target_information, fastq_dir, amplicons_dir)
    att_sequence_dirs = MEASURE_INTEGRATION(target_information, alignment_dir)
    GATHER_QC_INFO(trim_outputs[1], fastq_dir)
    cs2_attL_attR_dirs = RUN_CS2_FOR_INDELS(att_sequence_dirs[1], att_sequence_dirs[2], amplicons_dir, target_information)
    indel_tables = GET_INDEL_INFO_FROM_CS2_OUTPUT(cs2_attL_attR_dirs[0], cs2_attL_attR_dirs[1])
    COMBINE_INTEGRATION_AND_INDEL_INFO(att_sequence_dirs[0],indel_tables[0],indel_tables[1])
}
