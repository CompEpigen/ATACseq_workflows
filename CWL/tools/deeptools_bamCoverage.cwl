cwlVersion: v1.0
class: CommandLineTool
requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 20000
    #tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/deeptools:3.1.1

baseCommand: ["bamCoverage"]
arguments:  
  - valueFrom: $(inputs.bam.nameroot + ".bigwig")
    prefix: --outFileName
    position: 10
  - valueFrom: "bigwig"
    prefix: --outFileFormat
    position: 10
  - valueFrom: |
        ${ 
          if( inputs.spike_in_count == null ){
            return "RPGC"
          }
          else{
            return null 
          }
        }
    prefix: --normalizeUsing
    position: 10
  
inputs:
  bam:
    doc: bam file as input; needs bai index file in the same directory
    type: File
    secondaryFiles: .bai
    inputBinding:
        position: 100
        prefix: --bam
  effective_genome_size:
    doc: |
      the effectively mappable genome size, 
      see: https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
    type: long
    inputBinding:
        position: 10
        prefix: --effectiveGenomeSize
  bin_size:
    type: int
    default: 10
    inputBinding:
      prefix: --binSize
      position: 10
  ignoreForNormalization:
    doc: |
      List of space-delimited chromosome names that shall be ignored
      when calculating the scaling factor. 
    type: string?
    default: "chrX chrY chrM"
    inputBinding:
      prefix: --ignoreForNormalization
      position: 10
  spike_in_count:
    doc: number of reads aligned to the spike in genome, optional
    type: long?
    inputBinding:
      position: 10
      prefix: --scaleFactor
      valueFrom: |
        ${ 
          if( self == null ){
            return null
          }
          else{
            return (1.0 / parseFloat(self)) 
          }
        }

outputs:
  bigwig:
    type: File
    outputBinding:
      glob: $(inputs.bam.nameroot + ".bigwig")
    