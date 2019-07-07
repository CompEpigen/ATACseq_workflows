cwlVersion: v1.0
class: CommandLineTool
hints:
  ResourceRequirement:
    coresMin: 4
    ramMin: 15000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7
    
baseCommand: ["samtools", "flagstat"]
stdout: $(inputs.bam.nameroot + inputs.output_suffix)

inputs:
  bam:
    type: File
    inputBinding:
      position: 2
  output_suffix:
    type: string
    default: .flagStat
    
outputs:
  flagstat_output:
    type: stdout

