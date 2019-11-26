cwlVersion: v1.0
class: Workflow

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}
  MultipleInputFeatureRequirement: {}

inputs:
  sample_id:
    type: string
  fastq1:
    type: 
      type: array
      items: File
  fastq2: 
    type:
      type: array
      items: File
  reference:
    type: File
    secondaryFiles:
      - .fai
      - ^.1.bt2
      - ^.2.bt2
      - ^.3.bt2
      - ^.4.bt2
      - ^.rev.1.bt2
      - ^.rev.2.bt2
  reference_info:
    type: File
  adapter1: 
    type: string?
  adapter2:
    type: string?
  max_mapping_insert_length:
    type: long
    default: 2500
  macs2_genome_size:
    type: string
  macs2_qvalue:
    type: float
    default: 0.05
  effective_genome_size:
    type: long
  bin_size:
    type: int
    default: 10
  ignoreForNormalization:
    doc: |
      List of space-delimited chromosome names that shall be ignored
      when calculating the scaling factor. 
    type: string
    default: "chrX chrY chrM"
  
steps:
  trim_and_map:
    run: "../workflow_modules/trim_and_map.cwl"
    scatter: [fastq1, fastq2]
    scatterMethod: 'dotproduct'
    in:
      fastq1:
        source: fastq1
      fastq2: 
        source: fastq2
      reference:
        source: reference
      adapter1: 
        source: adapter1
      adapter2:
        source: adapter2
      is_paired_end:
        default: true
      max_mapping_insert_length:
        source: max_mapping_insert_length
    out:
      - raw_fastqc_zip
      - raw_fastqc_html
      - fastq1_trimmed
      - fastq2_trimmed
      - trim_galore_log
      - trimmed_fastqc_html
      - trimmed_fastqc_zip
      - bam
      - bowtie2_log
 
  merge_duprem_filter:
    run: "../workflow_modules/merge_duprem_filter.cwl"
    in:
      sample_id:
        source: sample_id
      bams: 
        source: trim_and_map/bam
      is_paired_end:
        default: true
    out:
      - duprem_fastqc_zip
      - duprem_fastqc_html
      - merged_flagstat_output
      - filtered_flagstat_output
      - duprem_flagstat_output
      - picard_markdup_log
      - bam

  name_sorting_filtered_bam:
      doc: samtools sort - sorting of filtered bam file by read name
      run: "../tools/samtools_sort_name.cwl"
      in:
        bam_unsorted:
          source: merge_duprem_filter/bam
      out:
        - bam_sorted
       
  converting_bam_to_bedpe:
    doc: bedtools bamtobed
    run: "../tools/bedtools_bamtobed_pe.cwl"
    in:
      bam:
        source: name_sorting_filtered_bam/bam_sorted
    out:
      - bedpe

  generating_atac_signal_tags:
    doc: 
    run: "../tools/generate_atac_signal_tags.cwl"
    in:
      bedpe_alignm: # paired alignments in bedpe format
        source: converting_bam_to_bedpe/bedpe
      output_basename:
        source: sample_id
    out:
      - bed_tn5_center_29bp
      - bed_tn5_center_73bp
      - bed_tn5_center_200bp
      - bed_tn5_center_1bp
      - bed_tn5_center_fragment
      - fragment_sizes_tsv
      - filtering_stats_tsv
      - frag_size_stats_tsv
      - irreg_mappings_bedpe

  generating_coverage_tracks:
    doc:
    run: "../workflow_modules/bed_to_coverage_track.cwl"
    scatter: [bed]
    scatterMethod: dotproduct
    in:
      bed:
        source: 
          - generating_atac_signal_tags/bed_tn5_center_29bp
          - generating_atac_signal_tags/bed_tn5_center_73bp
          - generating_atac_signal_tags/bed_tn5_center_200bp
          - generating_atac_signal_tags/bed_tn5_center_1bp
          - generating_atac_signal_tags/bed_tn5_center_fragment
      reference_info:
        source: reference_info
      effective_genome_size:
        source: effective_genome_size
      bin_size:
        source: bin_size
      ignoreForNormalization:
        source: ignoreForNormalization
    out:
      - bigwig
      - bam

  peak_calling_macs2_broad:
    doc: peak calling using macs2
    run: "../tools/macs2_callpeak_atac.cwl"
    scatter: [treatment_bed]
    scatterMethod: dotproduct
    in:
      treatment_bed:
        source:
          - generating_atac_signal_tags/bed_tn5_center_73bp
          - generating_atac_signal_tags/bed_tn5_center_200bp
          - generating_atac_signal_tags/bed_tn5_center_fragment
        linkMerge: merge_flattened
      genome_size:
        source: macs2_genome_size
      broad:
        valueFrom: ${return(true)}
      qvalue:
        source: macs2_qvalue
    out: 
      - peaks_bed
      - peaks_xls
      
  peak_calling_macs2_narrow:
    doc: peak calling using macs2
    run: "../tools/macs2_callpeak_atac.cwl"
    in:
      treatment_bed:
        source: generating_atac_signal_tags/bed_tn5_center_29bp
      genome_size:
        source: macs2_genome_size
      broad:
        valueFrom: ${return(false)}
      qvalue:
        source: macs2_qvalue
    out: 
      - peaks_bed
      - peaks_xls

  plot_fragment_size_distribution:
    run: "../tools/plot_frag_size_distr.cwl"
    in:
      fragment_sizes_tsv:
        source: generating_atac_signal_tags/fragment_sizes_tsv
      output_basename:
        source: sample_id
    out:
      - frag_size_distr_plot
      - frag_size_distr_tsv

  qc_plot_fingerprint:
    run: "../tools/deeptools_plotFingerprint.cwl"
    in:
      bam:
        source: merge_duprem_filter/bam
      sample_id:
        source: sample_id
      is_paired_end:
        default: true
    out:
      - qc_plot_fingerprint_plot  
      - qc_plot_fingerprint_tsv
      - qc_plot_fingerprint_stderr

  qc_phantompeakqualtools:
    run: "../tools/phantompeakqualtools.cwl"
    in:
      bam:
        source: merge_duprem_filter/bam
    out:
      - qc_crosscorr_summary  
      - qc_crosscorr_plot
      - qc_phantompeakqualtools_stderr
      - qc_phantompeakqualtools_stdout

  create_summary_qc_report:
    doc: |
      multiqc summarizes the qc results from fastqc 
      and other tools
    run: "../tools/multiqc_hack.cwl"
    in:
      qc_files_array_of_array:
        source:
          - trim_and_map/raw_fastqc_zip
          - trim_and_map/raw_fastqc_html
          - trim_and_map/trimmed_fastqc_html
          - trim_and_map/trimmed_fastqc_zip
          - trim_and_map/trim_galore_log
          - peak_calling_macs2_broad/peaks_bed
        linkMerge: merge_flattened
      qc_files_array:
        source:
          - trim_and_map/bowtie2_log
          - merge_duprem_filter/duprem_fastqc_zip
          - merge_duprem_filter/duprem_fastqc_html
          - generating_atac_signal_tags/frag_size_stats_tsv
          - generating_atac_signal_tags/fragment_sizes_tsv
          - generating_atac_signal_tags/filtering_stats_tsv
          - peak_calling_macs2_broad/peaks_xls
          - peak_calling_macs2_narrow/peaks_bed
          - peak_calling_macs2_narrow/peaks_xls
          - plot_fragment_size_distribution/frag_size_distr_tsv
          - qc_plot_fingerprint/qc_plot_fingerprint_tsv
          - qc_phantompeakqualtools/qc_phantompeakqualtools_stdout
          - qc_phantompeakqualtools/qc_crosscorr_summary
          - merge_duprem_filter/picard_markdup_log
        linkMerge: merge_flattened
      report_name:
        source: sample_id
    out:
      - multiqc_zip
      - multiqc_html

outputs:
  raw_fastqc_zip:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/raw_fastqc_zip
  raw_fastqc_html:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/raw_fastqc_html
  trim_galore_log:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/raw_fastqc_zip
  trimmed_fastqc_html:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/trimmed_fastqc_html
  trimmed_fastqc_zip:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: trim_and_map/trimmed_fastqc_zip
  bowtie2_log:
    type:
      type: array
      items: File
    outputSource: trim_and_map/bowtie2_log

  duprem_fastqc_zip:
    type:
      type: array
      items: File
    outputSource: merge_duprem_filter/duprem_fastqc_zip
  duprem_fastqc_html:
    type:
      type: array
      items: File
    outputSource: merge_duprem_filter/duprem_fastqc_html
  merged_flagstat_output:
    type: File
    outputSource: merge_duprem_filter/merged_flagstat_output
  filtered_flagstat_output:
    type: File
    outputSource: merge_duprem_filter/filtered_flagstat_output
  duprem_flagstat_output:
    type: File
    outputSource: merge_duprem_filter/duprem_flagstat_output
  bam:
    type: File
    secondaryFiles: .bai
    outputSource: merge_duprem_filter/bam
  picard_markdup_log:
    type: File
    outputSource: merge_duprem_filter/picard_markdup_log

  frag_size_stats_tsv:
    type: File
    outputSource: generating_atac_signal_tags/frag_size_stats_tsv
  filtering_stats_tsv:
    type: File
    outputSource: generating_atac_signal_tags/filtering_stats_tsv
  fragment_sizes_tsv:
    type: File
    outputSource: generating_atac_signal_tags/fragment_sizes_tsv
  irreg_mappings_bedpe:
    type: File
    outputSource: generating_atac_signal_tags/irreg_mappings_bedpe

  bam_signal_tags:
    type: 
      type: array
      items: File
    outputSource: generating_coverage_tracks/bam 
  bigwig_signal_tags:
    type: 
      type: array
      items: File
    outputSource: generating_coverage_tracks/bigwig 

  peaks_bed_macs2_broad:
    type:
      type: array
      items: 
        type: array
        items: File
    outputSource: peak_calling_macs2_broad/peaks_bed
  peaks_xls_macs2_broad:
    type:
      type: array
      items: File
    outputSource: peak_calling_macs2_broad/peaks_xls
  peaks_bed_macs2_narrow:
    type: 
      type: array
      items: File
    outputSource: peak_calling_macs2_narrow/peaks_bed
  peaks_xls_macs2_narrow:
    type: File
    outputSource: peak_calling_macs2_narrow/peaks_xls

  frag_size_distr_plot:
    type: File
    outputSource: plot_fragment_size_distribution/frag_size_distr_plot
  frag_size_distr_tsv:
    type: File
    outputSource: plot_fragment_size_distribution/frag_size_distr_tsv
  qc_plot_fingerprint_plot:
    type: File?
    outputSource: qc_plot_fingerprint/qc_plot_fingerprint_plot
  qc_plot_fingerprint_tsv:
    type: File?
    outputSource: qc_plot_fingerprint/qc_plot_fingerprint_tsv
  qc_plot_fingerprint_stderr:
    type: File
    outputSource: qc_plot_fingerprint/qc_plot_fingerprint_stderr
  qc_crosscorr_summary:
    type: File?
    outputSource: qc_phantompeakqualtools/qc_crosscorr_summary
  qc_crosscorr_plot:
    type: File?
    outputSource: qc_phantompeakqualtools/qc_crosscorr_plot
  qc_phantompeakqualtools_stderr:
    type: File?
    outputSource: qc_phantompeakqualtools/qc_phantompeakqualtools_stderr

  multiqc_zip:
    type: File
    outputSource: create_summary_qc_report/multiqc_zip
  multiqc_html:
    type: File
    outputSource: create_summary_qc_report/multiqc_html