cwlVersion: v1.0
class: CommandLineTool
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  InitialWorkDirRequirement:
    # This step is necessary since the input files
    # must be loaded into the working directory as there
    # is no way to specify the input file directly on the
    # command line.
    listing: |
      ${// script merges the to input arrays
        // into one array that fulfills the type 
        // requirement for "listing", which is
        // "{type: array, items: [File, Directory]}"

        var qc_files_array = inputs.qc_files_array;
        var qc_files_array_of_array = inputs.qc_files_array_of_array;
        var output_array = [];

        // add items of the qc_files_array to the output_array
        if ( qc_files_array != null ){
          for (var i=0; i<qc_files_array.length; i++){
            output_array.push(qc_files_array[i])
          }
        }

        // add items of the qc_files_array_of_array to the output_array
        if ( qc_files_array_of_array != null ){
          for (var i=0; i<qc_files_array_of_array.length; i++){ 
            for (var ii=0; ii<qc_files_array_of_array[i].length; ii++){
              output_array.push(qc_files_array_of_array[i][ii])
            }
          }
        }

        return output_array
      }

hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/multiqc:1.7
  

baseCommand: ["multiqc"]
arguments:
  - valueFrom: --zip-data-dir
    position: 1
  - valueFrom: "'log_filesize_limit: 100000000'"
    position: 1
    prefix: --cl_config
  - valueFrom: $(runtime.outdir)
    position: 2
    prefix: --outdir
  - valueFrom: $(runtime.outdir)
    position: 4
  
inputs:
  qc_files_array:
    doc: |
      qc files which shall be part of the multiqc summary;
      optional, only one of qc_files_array or qc_files_array_of_array 
      must be provided
    type:
      - "null"
      - type: array
        items: File
  qc_files_array_of_array:
    doc: |
      qc files which shall be part of the multiqc summary;
      optional, only one of qc_files_array or qc_files_array_of_array 
      must be provided
    type:
      - "null"
      - type: array
        items: 
          type: array
          items: File
  report_name:
    doc: name used for the html report and the corresponding zip file
    type: string
    default: multiqc
    inputBinding:
      valueFrom: $(self + "_report")
      prefix: --filename
      position: 3
      
outputs:
  multiqc_zip:
    type: File
    outputBinding:
      glob: $(inputs.report_name + "_report_data.zip")
  multiqc_html:
    type: File
    outputBinding:
      glob: $(inputs.report_name + "_report.html")
  