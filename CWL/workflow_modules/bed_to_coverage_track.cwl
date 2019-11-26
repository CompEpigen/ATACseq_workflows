cwlVersion: v1.0
class: Workflow

### INPUT PART:
##################################################
inputs:
  bed:
    type: File
  reference_info:
    type: File
  effective_genome_size:
    type: long
  bin_size:
    type: int
  ignoreForNormalization:
    doc: |
      List of space-delimited chromosome names that shall be ignored
      when calculating the scaling factor. 
    type: string
    default: "chrX chrY chrM"

### WORKFLOW STEPS:
##################################################
steps:
  converting_bed_to_bam:
    doc: |
      bedtools bedtobam - converts bed to bam;
      as most tools can handle bam as input but not always
      the bed format(e.g. deeptools), moreover, bam is compressed;
      therefore, it will be used as final output (instead of the bed file)
    run: "../tools/bedtools_bedtobam.cwl"
    in:
      bed:
        source: bed
      reference_info:
        source: reference_info
    out:
      - bam

  sorting_bam:
    doc: samtools sort - sorting of merged bam
    run: "../tools/samtools_sort.cwl"
    in:
      bam_unsorted:
        source: converting_bed_to_bam/bam
    out:
       - bam_sorted
       
  indexing_bam:
    doc: |
      samtools index - indexes sorted bam
    run: "../tools/samtools_index_hack.cwl"
    in:
      bam_sorted:
        source: sorting_bam/bam_sorted
    out:
       - bam_sorted_indexed

  converting_bam_to_bigwig:
    doc: |
      deeptools bamCoverage
    run: "../tools/deeptools_bamCoverage.cwl"
    in:
      bam:
        source: indexing_bam/bam_sorted_indexed
      effective_genome_size:
        source: effective_genome_size
      bin_size:
        source: bin_size
      ignoreForNormalization:
        source: ignoreForNormalization
    out:
      - bigwig


### OUTPUTS:
##################################################
outputs:
  bam:
    type: File
    outputSource: indexing_bam/bam_sorted_indexed
  bigwig:
    type: File
    outputSource: converting_bam_to_bigwig/bigwig
    
  
    
    
    