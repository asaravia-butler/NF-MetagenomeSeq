#!/usr/bin/env nextflow
nextflow.enable.dsl = 2




// process to indentify antibiotic genes
process RGI {

    tag "Running rgi on ${sample_id}-s reads..." 
    label "read_based_outputs"
    input:
        tuple val(sample_id), path(reads), val(isPaired)
        each path(card_json) 
    output:
        tuple val(sample_id), path("${sample_id}.overall_mapping_stats.txt"), emit: stats
        tuple val(sample_id), path("${sample_id}_mapping_summary.csv") , emit: summary
        path("versions.txt"), emit: version


    script:
      """

      rgi card_annotation -i ${card_json}

      rgi load --card_json ${card_json} \\
              --debug --local --card_annotation card_database_v4.0.0.fasta \\
              --card_annotation_all_models card_database_v4.0.0_all.fasta

      if [ ${isPaired} == true ];then

        rgi bwt --read_one ${reads[0]} \\
            --read_two ${reads[1]} \\
            --output_file ${sample_id} \\
            --local \\
            --aligner bowtie2 \\
            --threads ${task.cpus}
      else

        rgi bwt --read_one ${reads[0]} \\
            --output_file ${sample_id}_rgi_bwt \\
            --local \\
            --aligner bowtie2 \\
            --threads ${task.cpus}

      fi


      percent=`grep 'Mapped reads' ${sample_id}.overall_mapping_stats.txt | grep -oP '\\(\\K[^\\)]+' | tr -d '%'`
      echo ${sample_id},\${percent} > ${sample_id}_mapping_summary.csv
      echo "RGI 6.0.3" > versions.txt
      """

}
