cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: genomicpariscentre/macs2:2.1.0.20140616
 
### BASE COMMAND AND ARGUMENTS:
##################################################
baseCommand: ["macs2", "callpeak"]
arguments:  
  - valueFrom: "BED"
    prefix: "--format"
    position: 1
  - valueFrom: "--nomodel"
    position: 2
  - valueFrom: "all"
    prefix: "--keep-dup"
    position: 2
  - valueFrom: $(inputs.treatment_bed.nameroot + ".macs2")
    prefix: "--name"
    position: 100

inputs:
  treatment_bed:
    type: File
    inputBinding:
        position: 101
        prefix: "--treatment"
  genome_size:
    doc: can be "mm", "hs", "ce", "dm", or the total number of genomic bp 
    type: string
    inputBinding:
        position: 3
        prefix: "--gsize"
  broad:
    type: boolean
    inputBinding:
        position: 3
        prefix: "--broad"
  qvalue:
    type: float
    inputBinding:
        position: 3
        prefix: "--qvalue"

 
outputs:
  peaks_bed:    
    type: 
      type: array
      items: File
    outputBinding:
      glob: "*Peak"
  peaks_xls:
    type: File
    outputBinding:
      glob: "*_peaks.xls"