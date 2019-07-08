cwlVersion: v1.0
class: CommandLineTool
hints:
  ResourceRequirement:
    coresMin: 1
    ramMin: 15000
    #tmpdirMin: 10000
  DockerRequirement:
    dockerPull: kerstenbreuer/samtools:1.7
requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  InitialWorkDirRequirement:
    listing: 
      - entryname: generate_atac_signal_tag.sh
        entry: | 
          BEDPE="\$1"
          OUTPUT_BASENAME="\$2"
          touch irreg_mappings.bedpe 
          touch fragment_sizes.txt
          touch tn5_center_29bp_unsorted.bed 
          touch tn5_center_73bp_unsorted.bed
          touch tn5_center_200bp_unsorted.bed 
          touch tn5_center_fragment_unsorted.bed
          touch tn5_center_1bp_unsorted.bed 
          awk generate_atac_signal_tag.awk "\${BEDPE}"
          LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n tn5_center_29bp_unsorted.bed > "\${OUTPUT_BASENAME}_tn5_center_29bp.bed"
          LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n tn5_center_73bp_unsorted.bed > "\${OUTPUT_BASENAME}_tn5_center_73bp.bed"
          LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n tn5_center_200bp_unsorted.bed > "\${OUTPUT_BASENAME}_tn5_center_200bp.bed"
          LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n tn5_center_fragment_unsorted.bed > "\${OUTPUT_BASENAME}_tn5_center_fragment.bed"
          LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n -k6,6n tn5_center_1bp_unsorted.bed > "\${OUTPUT_BASENAME}_tn5_center_1bp.bed"
          rm tn5_center_29bp_unsorted.bed 
          rm tn5_center_73bp_unsorted.bed   
          rm tn5_center_200bp_unsorted.bed  
          rm tn5_center_fragment_unsorted.bed 
          rm tn5_center_1bp_unsorted.bed 
          mv irreg_mappings.bedpe "\${OUTPUT_BASENAME}_irreg_mappings.bedpe"
          mv fragment_sizes.txt "\${OUTPUT_BASENAME}_fragment_sizes.txt"
      - entryname: generate_atac_signal_tag.awk
        entry: |
          function tn5_center_ext(seqname, name, score, bp, center1, center2, output_file) { 
            # bp must be uneven
            tag_start = (center1-(bp/2));
            tag_end = (center1+(bp/2)+1);
            if(tag_start < 0){tag_start=0};
            if(tag_end < 0){tag_end=0};
            print seqname, tag_start, tag_end, name, score > output_file;
            tag_start = (center2-(bp/2));
            tag_end = (center2+(bp/2)+1);
            if(tag_start < 0){tag_start=0};
            if(tag_end < 0){tag_end=0};
            print seqname, tag_start, tag_end, name"_mate", score > output_file;
          }
          BEGIN {
              OFS="\t";
              chrM_read_count=0;
              interchrom_map_read_count=0;
              regular_read_count=0;
              irregular_read_count=0;
              too_small_fragment_count=0;
              nucl_free_fragment_count=0;
              nucl_bound_fragment_count=0;
              wrong_strand_orient_count=0;
          }
          {
            if ( $1=="chrM" || $4=="chrM") {
                  irregular_read_count += 2;
                  if ( $1==$4 ) {
                      chrM_read_count += 2;
                      print $0, "chrM" > "irreg_mappings.bedpe" ;
                  }
                  else {
                      chrM_read_count += 1;
                      interchrom_map_read_count += 2;
                      print $0, "interchrom_map" > "irreg_mappings.bedpe";
                  }
              }
              else if ( $1!=$4 ) {
                  irregular_read_count += 2;
                  interchrom_map_read_count += 2;
                  print $0, "interchrom_map" > "irreg_mappings.bedpe";
              }
              else if ( $9==$10 ) {
                  irregular_read_count += 2;
                  wrong_strand_orient_count += 2;
                  print $0, "wrong_read_pair_orientation" > "irreg_mappings.bedpe";
              }
              else {
                  right_orientation=0;
                  # shift to center and calc fragment size
                  if ( $9=="+" && $2<=$5 && $3<=$6 ){
                      right_orientation=1;
                      fragment_size=$6-$2;
                      # first read:
                      center1=$2+4;
                      # second read
                      center2=$6-5;
                  }
                  else if ( $9=="-" && $2>=$5 && $3>=$6 ){
                      right_orientation=1;
                      fragment_size=$3-$5;
                      # first read:
                      center1=$5+4;
                      # second read
                      center2=$3-5;
                  }
                  else {
                      wrong_strand_orient_count += 2;
                      print $0, "wrong_read_pair_orientation" > "irreg_mappings.bedpe";
                  }
                  print "fragment_size" >  "fragment_sizes.txt";
                  if ( fragment_size>=38 && right_orientation ) {
                      regular_read_count += 2;
                      tn5_center_ext($1, $7, $8, 29, center1, center2, "tn5_center_29bp_unsorted.bed");
                      tn5_center_ext($1, $7, $8, 73, center1, center2, "tn5_center_73bp_unsorted.bed");
                      tn5_center_ext($1, $7, $8, 200, center1, center2, "tn5_center_200bp_unsorted.bed");
                      tn5_center_ext($1, $7, $8, 1, center1, center2, "tn5_center_1bp_unsorted.bed");
                      print $1, center1, center2, $7, $8 >  "tn5_center_fragment_unsorted.bed";
                      if (fragment_size>=185) { 
                          nucl_bound_fragment_count += 1;
                      }
                      else {
                          nucl_free_fragment_count += 1;
                      }
                  }
                  else {
                      if ( right_orientation ){
                          irregular_read_count += 2;
                          too_small_fragment_count += 1;
                          print $0, "too_small_frag_size" > "irreg_mappings.bedpe" ;
                      }
                  }
            }
          }
          END {
              print "\# id: \"Filtering Statistics\"" >  "filtering_stats_mqc.tsv";
              print "\# description: \"- This section shows statistics on read filtering: (1) reads pairs that are mapping to different chromosomes as well as (2) read pairs located on ChrM are filtered out; (3) read pairs that have a wrong orientation towards each other (e.g. both reads on same strand, or reads pointing to different direction) are removed, too. (4) Only the remaining regular reads are used for the fragment size analysis and the generation of atac signal tracks. In addition to the filtering shown here, reads were also selected for high mapping quality and to be mapped in a proper pair.\"" >  filtering_stats_mqc.tsv;
              print "\# plot_type: \"bargraph\"" >  "filtering_stats_mqc.tsv";
              print "chrM_reads", "chrM_read_count" >  "filtering_stats_mqc.tsv";
              print "interchrom_map_reads", "interchrom_map_read_count" >  "filtering_stats_mqc.tsv";
              print "wrong_read_pair_orientation", "wrong_strand_orient_count" >  "filtering_stats_mqc.tsv";
              print "regular_reads", "regular_read_count" >  "filtering_stats_mqc.tsv";
              print "\# id: \"Fragment Length Classification\"" >  "frag_size_classification_mqc.tsv";
              print "\# description: \"Fragments are classified by their size: (1) nucleosome free fragements are smaller than a typical a nucleosome binding region while (2)) potentially nucleosome bound fragments are larger;(3) fragments that are classified as too small are shorter than expected for a Tn5 digestion. For these calculations, the DNA span that is covered by the Tn5 enzyme during the trasposition is taken into account.\"" >  frag_size_classification_mqc.tsv;
              print "\# plot_type: \"bargraph\"" >  "frag_size_classification_mqc.tsv";
              print "too_small_fragments", "too_small_fragment_count" >  "frag_size_classification_mqc.tsv";
              print "nucl_free_fragments", "nucl_free_fragment_count" >  "frag_size_classification_mqc.tsv";
              print "pot_nucl_bound_fragments", "nucl_bound_fragment_count" >  "frag_size_classification_mqc.tsv";
          }


  
  
### BASE COMMAND AND ARGUMENTS:
##################################################
baseCommand: ["bash", "generate_atac_signal_tag.sh"]
        
### INPUT PART:
##################################################
inputs:
  bedpe_alignm:
    type: File
    inputBinding:
      position: 1
  output_basename:
    type: string
    inputBinding:
      position: 2
 
### OUTPUT PART:
##################################################
outputs:
  bed_tn5_center_29bp:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_tn5_center_29bp.bed")
  bed_tn5_center_73bp:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_tn5_center_73bp.bed")
  bed_tn5_center_200bp:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_tn5_center_200bp.bed")
  bed_tn5_center_1bp:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_tn5_center_1bp.bed")
  bed_tn5_center_fragment:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_tn5_center_fragment.bed")
  # bed_signal_tags:
  #   type:
  #     type: array
  #     items: File
  #   outputBinding:
  #     glob: "*.bed"
    
  fragment_sizes_tsv:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_fragment_sizes.txt")
  filtering_stats_tsv:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_filtering_stats_mqc.tsv")
  frag_size_stats_tsv:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_frag_size_classification_mqc.tsv")
  irreg_mappings_bedpe:
    type: File
    outputBinding:
      glob: $(inputs.output_basename + "_irreg_mappings.bedpe")