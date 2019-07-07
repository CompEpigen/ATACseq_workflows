
cwlVersion: v1.0
class: Workflow

requirements:
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
  InlineJavascriptRequirement: {}

### INPUT PART:
##################################################
inputs:
  sample_id:
    type: string
  bams:
    type:
      type: array
      items: File
  is_paired_end:
    type: boolean
        
### WORKFLOW STEPS:
##################################################
steps:
  lane_replicate_merging:
    doc: samtools merge - merging bam files of lane replicates
    run: "../tools/samtools_merge.cwl"
    in:
      bams:
        source: bams
      output_name:
        source: sample_id
        valueFrom: $(self + ".bam")
    out:
       - bam_merged

  sorting_merged_bam:
    doc: samtools sort - sorting of merged bam
    run: "../tools/samtools_sort.cwl"
    in:
      bam_unsorted:
        source: lane_replicate_merging/bam_merged
    out:
       - bam_sorted
       
  flagstat_merged:
    doc: samtools flagstat on merged bams
    run: "../tools/samtools_flagstat.cwl"
    in:
      bam: sorting_merged_bam/bam_sorted
    out:
      - flagstat_output

  filter_by_mapq:
    doc: samtools view
    run: "../tools/samtools_view_filter.cwl"
    in:
      bam:
        source: sorting_merged_bam/bam_sorted
      is_paired_end:
        source: is_paired_end
    out:
      - bam_filtered

  flagstat_filtered:
    doc: samtools flagstat on quality filtered bams
    run: "../tools/samtools_flagstat.cwl"
    in:
      bam: filter_by_mapq/bam_filtered
    out:
      - flagstat_output

  remove_duplicates:
    doc: picard markdup - emoves duplicates from a single sorted bam file.
    run: "../tools/picard_markdup.cwl"
    in:
      bam_sorted:
        source: filter_by_mapq/bam_filtered
    out:
      - bam_duprem
      - picard_markdup_log

  flagstat_duprem:
    doc: samtools flagstat on merged bams
    run: "../tools/samtools_flagstat.cwl"
    in:
      bam: remove_duplicates/bam_duprem
    out:
      - flagstat_output

  indexing_duprem_bam:
    doc: |
      samtools index - indexes sorted bam
    run: "../tools/samtools_index_hack.cwl"
    in:
      bam_sorted:
        source: remove_duplicates/bam_duprem
    out:
       - bam_sorted_indexed

  qc_duprem:
    doc: fastqc - quality control for reads directly after mapping
    run: "../tools/fastqc.cwl"
    in:
      bam:
        source: indexing_duprem_bam/bam_sorted_indexed
    out:
      - fastqc_zip
      - fastqc_html
      
### OUTPUTS:
##################################################
outputs:
  #bam_merged:
  #  type: File
  #  outputSource: lane_replicate_merging/bam_merged
  #bam_merged_sorted:
  #  type: File
  #  outputSource: sorting_merged_bam/bam_sorted
  #bam_merged_duprem:
  #  type: File
  #  outputSource: remove_duplicates/bam_duprem
  #bam_merged_duprem_filtered:
  #  type: File
  #  outputSource: filter_by_mapq/bam_filtered
  duprem_fastqc_zip:
    type: 
      type: array
      items: File
    outputSource: qc_duprem/fastqc_zip
  duprem_fastqc_html:
    type: 
      type: array
      items: File
    outputSource: qc_duprem/fastqc_html
  merged_flagstat_output:
    type: File
    outputSource: flagstat_merged/flagstat_output
  filtered_flagstat_output:
    type: File
    outputSource: flagstat_filtered/flagstat_output
  duprem_flagstat_output:
    type: File
    outputSource: flagstat_duprem/flagstat_output
  bam:
    type: File
    secondaryFiles: .bai
    outputSource: indexing_duprem_bam/bam_sorted_indexed
  picard_markdup_log:
    type: File
    outputSource: remove_duplicates/picard_markdup_log
    
  
    
    
    