nextflow.enable.dsl=2

// params.r1 = "/data/reads/CM_Custom_TB000175c_20230613/4L-550-rep1_L1_ds.228ab3596da04f878fcc010e9dc007b5/4L-550-rep1_S1_L001_R1_001.fastq.gz"
// params.r2 = "/data/reads/CM_Custom_TB000175c_20230613/4L-550-rep1_L1_ds.228ab3596da04f878fcc010e9dc007b5/4L-550-rep1_S1_L001_R2_001.fastq.gz"
// params.metadata = "/data/CM_Custom_TB000175c_20230613_analysis/20230329-final_panel_up.csv"
params.reference = "/data/references/hg38.fa"
params.outdir = "/data/4L-550-rep1/"
// params.sample_name = "test"
params.ATTP_REG = 'GTGGTTTGTCTGGTCAACCACCGCGGT'
params.ATTP_PRIME = 'CTCAGTGGTGTACGGTACAAACCCA'
params.deduplication_method = 'STRICT'
params.collapse_condition = 'Complete'
params.project_config_file = ''

process ADAPTER_AND_POLY_G_TRIM {
    maxForks 1
    publishDir "${params.outdir}/${sample_name}/qc/", pattern: 'fastp.json'

    input:
        tuple val(sample_name), path(R1), path(R2), path(metadata)
    output:
        path "trimmed.fastq.bgz"
        path "fastp.json"

    script:
        """

        fastp -m -c --include_unmerged --overlap_len_require 10 -i ${R1} -I ${R2} --merged_out trimmed.fastq.gz -w 16 -g
        gunzip -c trimmed.fastq.gz | bgzip -@ 8 > trimmed.fastq.bgz
        """
}

process ALIGN_READS {
    maxForks 1
    publishDir "${params.outdir}/${sample_name}/"

    input:
        tuple val(sample_name), path(R1), path(R2), path(metadata)
        path reference
        path fastq
        val dedup_method

    output:
        tuple path("initial_alignment.bam"), path("initial_alignment.bam.bai")
        path "${reference}.fai"

    script:
    if (dedup_method == "STRICT")
        """
        minimap2 -ax sr -t 16 ${reference} ${fastq} > initial_alignment.sam 
        samtools view -@ 16 -b initial_alignment.sam > initial_alignment.bam
        samtools sort -@ 16 initial_alignment.bam > initial_alignment_sorted.bam
        mv initial_alignment_sorted.bam initial_alignment.bam
        rm initial_alignment.sam
        samtools sort -@ 16 -n -o namesort.bam initial_alignment.bam
        samtools fixmate -m namesort.bam fixmate.bam
        samtools sort -@ 16  -o positionsort.bam fixmate.bam
        samtools markdup --barcode-rgx ":(\\w+)\$"  -r positionsort.bam markdup.bam  
        mv markdup.bam initial_alignment.bam
        samtools index initial_alignment.bam
        samtools faidx ${reference}
        rm positionsort.bam
        rm fixmate.bam
        rm namesort.bam
        """
    else if (dedup_method == "LOOSE")
        """
        minimap2 -ax sr -t 16 ${reference} ${fastq} > initial_alignment.sam 
        samtools view -@ 16 -b initial_alignment.sam > initial_alignment.bam
        samtools sort -@ 16 initial_alignment.bam > initial_alignment_sorted.bam
        mv initial_alignment_sorted.bam initial_alignment.bam
        rm initial_alignment.sam
        samtools sort -@ 16 -n -o namesort.bam initial_alignment.bam
        samtools fixmate -m namesort.bam fixmate.bam
        samtools sort -@ 16  -o positionsort.bam fixmate.bam
        samtools markdup  -r positionsort.bam markdup.bam
        mv markdup.bam initial_alignment.bam
        samtools index initial_alignment.bam
        samtools faidx ${reference}
        rm positionsort.bam
        rm fixmate.bam
        rm namesort.bam
        """
    }

process GET_TARGET_INFORMATION {
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
    maxForks 1
    publishDir "${params.outdir}/${sample_name}/"
    input:
        tuple val(sample_name), path(R1), path(R2), path(metadata_fn)
        path target_info
        tuple path(bam_file), path(bam_file_index)
        path fastq
    output:
        path "*extracted_reads"
        path "*read_counts_per_site.csv"
    
    script:
    """
    extract_target_reads.py --target_info ${target_info} --bam_file ${bam_file} --sample_name ${sample_name}
    """
}

process GENERATE_AMPLICONS {
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
    maxForks 1
    publishDir "${params.outdir}/${sample_name}/"
    input:
        tuple val(sample_name), path(R1), path(R2), path(metadata)
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
    publishDir "${params.outdir}/${sample_name}", pattern: 'integration_stats.csv'
    input:
        tuple val(sample_name), path(R1), path(R2), path(metadata)
        path target_info
        path alignment_dir
    output:
        path "*integration_stats.csv" 
        path "*attL_extracted_reads" 
        path "*attR_extracted_reads" 
        path "*beacon_extracted_reads" 
    
    script:
    """
    compute_integration_percentage.py --target_info ${target_info} --alignment_dir ${alignment_dir} --sample_name ${sample_name}
    """
}

process GATHER_QC_INFO {
    publishDir "${params.outdir}/${sample_name}/qc"
    input:
        tuple val(sample_name), path(R1), path(R2), path(metadata)
        path json_file
        path fastq_dir
    output:
        path '*qc_summary.csv'
    script:
    """
    gather_qc_stats.py --json_file ${json_file} --fastq_dir ${fastq_dir} --sample_name ${sample_name}
    """
}

process RUN_CS2_FOR_INDELS {
    maxForks 1
    publishDir "${params.outdir}/${sample_name}/cs2_output"
    input:
        tuple val(sample_name), path(R1), path(R2), path(metadata)
        path attL_sequences
        path attR_sequences
        path beacon_sequences
        path amplicon_dir
        path target_info
    output:
        path 'cs2_attL'
        path 'cs2_attR'
        path 'cs2_beacon'
   
    script:
    """
    run_cs2.py --amplicon_dir ${amplicon_dir} --target_info ${target_info}
    """
}

process GET_INDEL_INFO_FROM_CS2_OUTPUT {
    maxForks 1
    publishDir "${params.outdir}/${sample_name}/indel_info"
    input:
        tuple val(sample_name), path(R1), path(R2), path(metadata)
        path cs2_attL
        path cs2_attR
    output:
        path "*attL_indel_table.csv"
        path "*attR_indel_table.csv"
   
    script:
    """
    get_indel_info.py --cs2_directory ${cs2_attL} --output_fn ${sample_name}_attL_indel_table.csv
    get_indel_info.py --cs2_directory ${cs2_attR} --output_fn ${sample_name}_attR_indel_table.csv
    """
}

process COMBINE_INTEGRATION_AND_INDEL_INFO {
    publishDir "${params.outdir}/${sample_name}"
    input:
        tuple val(sample_name), path(R1), path(R2), path(metadata)
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

process GENERATE_REPORT {  
  input:
    path project_config_file
    path integration_stats_files
    path attL_indel_table_files
    path attR_indel_table_files
    path read_counts_per_site_files
    path qc_summary_files
    path extracted_reads_dirs
    val collapse_condition
  output:
    path "*.xlsx"

  script:
    """
    collate_results.py --project_config_file ${project_config_file} --integration_stats_files ${integration_stats_files} --attL_indel_table_files ${attL_indel_table_files} --attR_indel_table_files ${attR_indel_table_files} --read_counts_per_site_files ${read_counts_per_site_files} --qc_summary_files ${qc_summary_files} --extracted_reads_dirs ${extracted_reads_dirs} --collapse_condition ${collapse_condition}
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

    Channel.fromPath(params.project_config_file)
    .splitCsv(header: true, sep: ',')
    .map { row -> 
        def r1 = "${row.fastq_dir}/*R1*.fastq.gz"
        def r2 = "${row.fastq_dir}/*R2*.fastq.gz"
        [
            row.sample_name, 
            file(r1), 
            file(r2), 
            file(row.probe_list)
        ] 
    }
    .set { input_ch }

    trimmed_and_merged_fastq = ADAPTER_AND_POLY_G_TRIM(input_ch)

    initial_alignment = ALIGN_READS(input_ch, params.reference, trimmed_and_merged_fastq[0], params.deduplication_method)

    probe_information = GET_TARGET_INFORMATION(input_ch, params.ATTP_REG, params.ATTP_PRIME)

    fastq_dir = EXTRACT_TARGET_READS(input_ch, probe_information, initial_alignment[0], trimmed_and_merged_fastq[0])

    amplicons_dir = GENERATE_AMPLICONS(input_ch, probe_information)

    alignment_dir = ALIGN_TARGET_READS(input_ch, probe_information, fastq_dir[0], amplicons_dir)

    att_sequence_dirs = MEASURE_INTEGRATION(input_ch, probe_information, alignment_dir)

    qc_summary = GATHER_QC_INFO(input_ch, trimmed_and_merged_fastq[1], fastq_dir[0])

    cs2_attL_attR_dirs = RUN_CS2_FOR_INDELS(input_ch, att_sequence_dirs[1], att_sequence_dirs[2], att_sequence_dirs[3], amplicons_dir, probe_information)

    indel_tables = GET_INDEL_INFO_FROM_CS2_OUTPUT(input_ch, cs2_attL_attR_dirs[0], cs2_attL_attR_dirs[1])

    COMBINE_INTEGRATION_AND_INDEL_INFO(input_ch, att_sequence_dirs[0], indel_tables[0], indel_tables[1])

    att_sequence_dirs[0].collect().view()
    
    GENERATE_REPORT(params.project_config_file,att_sequence_dirs[0].collect(), indel_tables[0].collect(), indel_tables[1].collect(), fastq_dir[1].collect(), qc_summary.collect(), fastq_dir[0].collect(), params.collapse_condition)
}
