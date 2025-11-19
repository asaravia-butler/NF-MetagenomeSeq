#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
//params.additional_filename_prefix = ""
//params.assay_suffix = "_GLmetagenomics"

/**************************************************************************************** 
********************* Read-based processing using Humann3 *******************************
****************************************************************************************/
// Make Kraken, Kaiju and Humann3 databases
include { SETUP_KAIJU; SETUP_KRAKEN; make_humann_db } from "./database_creation.nf"

include { METAPHLAN2KRONA; KRONA_REPORT as METAPHLAN_REPORT } from "./visualize_taxonomy.nf"
include { KRAKEN_CLASSIFY; KRAKEN2TABLE } from "./assign_taxonomy.nf"
include { KAIJU_CLASSIFY; KAIJU2TABLE } from "./assign_taxonomy.nf"
include { KRONA_REPORT as KRAKEN_REPORT } from "./visualize_taxonomy.nf"
include { KRAKEN2KRONA; KAIJU2KRONA } from "./visualize_taxonomy.nf"

include { KRONA_REPORT as KAIJU_REPORT } from "./visualize_taxonomy.nf"

/*
    This process runs humann3 and metaphlan4 on each individual sample generating the
    read-based functional annotations and taxonomic classifications.
*/

process HUMANN {

    tag "Running humann on ${sample_id}-s reads..."
    label "read_based"
    

    input:
        tuple val(sample_id), path(reads), val(isPaired)
        path(chocophlan_dir)
        path(uniref_dir)
        path(metaphlan_dir)
        
    output:
        path("${sample_id}-humann3-out-dir/${sample_id}_genefamilies.tsv"), emit: genefamilies
        path("${sample_id}-humann3-out-dir/${sample_id}_pathabundance.tsv"), emit: pathabundance
        path("${sample_id}-humann3-out-dir/${sample_id}_pathcoverage.tsv"), emit: pathcoverage
        tuple val(sample_id), path("${sample_id}-humann3-out-dir/${sample_id}_metaphlan_bugs_list.tsv"), emit: metaphlan_bugs_list
        path("versions.txt"), emit: version
    script:
        """
        zcat ${reads} > ${sample_id}-reads.tmp.fq

        humann --input ${sample_id}-reads.tmp.fq \\
                    --output ${sample_id}-humann3-out-dir/ \\
                    --threads ${task.cpus} \\
                    --output-basename ${sample_id} \\
                    --metaphlan-options "--bowtie2db ${metaphlan_dir} --unclassified_estimation --add_viruses --sample_id ${sample_id}" \\
                    --nucleotide-database ${chocophlan_dir} \\
                    --protein-database ${uniref_dir} \\
                    --bowtie-options "--sensitive --mm" && \\
        mv ${sample_id}-humann3-out-dir/${sample_id}_humann_temp/${sample_id}_metaphlan_bugs_list.tsv \\
               ${sample_id}-humann3-out-dir/${sample_id}_metaphlan_bugs_list.tsv

        humann3 --version  > versions.txt
        """
}


/*
    This process combines the read-based humann3 output functional 
    tables from indiviual samples into single tables across the GLDS dataset.
*/

process COMBINE_READ_BASED_PROCESSING_TABLES {

    tag "Combining the read based processing tables..."
    label "read_based"

    input:
        path(gene_families)
        path(path_abundances)
        path(path_coverages)
        path(utilities_path)
    output:
        path("${params.additional_filename_prefix}gene-families-initial.tsv"), emit: gene_families 
        path("${params.additional_filename_prefix}pathway-abundances-initial.tsv"), emit: path_abundances
        path("${params.additional_filename_prefix}pathway-coverages-initial.tsv"), emit: path_coverages
        path("versions.txt"), emit: version
    script:
        """
        if [ ${params.use_conda} == true ]; then
            # Setting humann3 utilities location (can be off if we pointed to
            # a previously installed database, and doesn't hurt to reset if it was already good-to-go)
            humann_config --update database_folders utility_mapping ${utilities_path} > /dev/null 2>&1
        fi

        # they each need to be in the same directories to be merged
        mkdir -p gene-family-results/ path-abundance-results/ path-coverage-results/
        cp ${gene_families} gene-family-results/ 
        cp ${path_abundances} path-abundance-results/
        cp ${path_coverages} path-coverage-results/

        humann_join_tables -i gene-family-results/ -o ${params.additional_filename_prefix}gene-families-initial.tsv > /dev/null 2>&1
        humann_join_tables -i path-abundance-results/ -o ${params.additional_filename_prefix}pathway-abundances-initial.tsv > /dev/null 2>&1
        humann_join_tables -i path-coverage-results/ -o ${params.additional_filename_prefix}pathway-coverages-initial.tsv > /dev/null 2>&1

        humann3 --version  > versions.txt
        """
}


/*
    The read-based functional annotation tables have taxonomic info and non-taxonomic info mixed
    together initially. Humann comes with utility scripts to split these. This process does that,
    generating non-taxonomically grouped functional info files and taxonomically grouped ones.
*/

process SPLIT_READ_BASED_PROCESSING_TABLES {

    tag "Splitting humann stratified tables..."
    label "read_based"
    label "read_based_outputs"

    input:
        path(gene_families)
        path(path_abundances)
        path(path_coverages)
    output:
        path("${params.additional_filename_prefix}Gene-families${params.assay_suffix}.tsv"), emit: gene_families
        path("${params.additional_filename_prefix}Gene-families-grouped-by-taxa${params.assay_suffix}.tsv"), emit: gene_families_grouped
        path("${params.additional_filename_prefix}Pathway-abundances${params.assay_suffix}.tsv"), emit: path_abundances
        path("${params.additional_filename_prefix}Pathway-abundances-grouped-by-taxa${params.assay_suffix}.tsv"), emit: path_abundances_grouped
        path("${params.additional_filename_prefix}Pathway-coverages${params.assay_suffix}.tsv"), emit: path_coverages 
        path("${params.additional_filename_prefix}Pathway-coverages-grouped-by-taxa${params.assay_suffix}.tsv"), emit: path_coverages_grouped
        path("versions.txt"), emit: version
    script:
        """
        [ -d temp_processing/ ] && rm -rf temp_processing/
        mkdir temp_processing/

        # Gene Families
        humann_split_stratified_table -i ${gene_families} -o temp_processing/ > /dev/null 2>&1
        mv temp_processing/${params.additional_filename_prefix}gene-families-initial_stratified.tsv \\
           ${params.additional_filename_prefix}Gene-families-grouped-by-taxa${params.assay_suffix}.tsv
        
        mv temp_processing/${params.additional_filename_prefix}gene-families-initial_unstratified.tsv \\
           ${params.additional_filename_prefix}Gene-families${params.assay_suffix}.tsv

        # Pathway Abundance
        humann_split_stratified_table -i ${path_abundances} -o temp_processing/ > /dev/null 2>&1
        mv temp_processing/${params.additional_filename_prefix}pathway-abundances-initial_stratified.tsv \\
           ${params.additional_filename_prefix}Pathway-abundances-grouped-by-taxa${params.assay_suffix}.tsv

        mv temp_processing/${params.additional_filename_prefix}pathway-abundances-initial_unstratified.tsv \\
           ${params.additional_filename_prefix}Pathway-abundances${params.assay_suffix}.tsv

        # Pathway Coverage
        humann_split_stratified_table -i ${path_coverages} -o temp_processing/ > /dev/null 2>&1
        mv temp_processing/${params.additional_filename_prefix}pathway-coverages-initial_stratified.tsv \\
           ${params.additional_filename_prefix}Pathway-coverages-grouped-by-taxa${params.assay_suffix}.tsv

        mv temp_processing/${params.additional_filename_prefix}pathway-coverages-initial_unstratified.tsv \\
           ${params.additional_filename_prefix}Pathway-coverages${params.assay_suffix}.tsv

        humann3 --version  > versions.txt
        """
}


/*
    This process generates some normalized tables of the read-based functional outputs from
    humann that are more readily suitable for across sample comparisons.
*/

process GEN_NORMALIZED_READ_BASED_PROCESSING_TABLES {

    tag "Generating normalized humann tables..."
    label "read_based"
    label "read_based_outputs"

    input:
       path(gene_families)
       path(path_abundances)
    output:
        path("${params.additional_filename_prefix}Gene-families-cpm${params.assay_suffix}.tsv"), emit: gene_families
        path("${params.additional_filename_prefix}Pathway-abundances-cpm${params.assay_suffix}.tsv"), emit: path_abundances
        path("versions.txt"), emit: version
    script:
        """
        humann_renorm_table \\
               -i ${gene_families} \\
               -o ${params.additional_filename_prefix}Gene-families-cpm${params.assay_suffix}.tsv \\
               --update-snames > /dev/null 2>&1

        humann_renorm_table \\
               -i ${path_abundances} \\
               -o ${params.additional_filename_prefix}Pathway-abundances-cpm${params.assay_suffix}.tsv \\
               --update-snames > /dev/null 2>&1

        humann3 --version  > versions.txt
        """
}


/*
    This process summarizes the read-based humann annotations based on Kegg Orthlogy terms.
*/

process GEN_READ_BASED_PROCESSING_KO_TABLE {

    tag "Retrieving Kegg Orthologs..."
    label "read_based"
    label "read_based_outputs"
    
    input:
        path(gene_families)
    output:
        path("${params.additional_filename_prefix}Gene-families-KO-cpm${params.assay_suffix}.tsv"), emit: gene_families
        path("versions.txt"), emit: version
    script:
        """
        humann_regroup_table \\
              -i ${gene_families} \\
              -g uniref90_ko 2> /dev/null | \\
        humann_rename_table \\
               -n kegg-orthology 2> /dev/null | \\
        humann_renorm_table \\
               -o ${params.additional_filename_prefix}Gene-families-KO-cpm${params.assay_suffix}.tsv \\
               --update-snames > /dev/null 2>&1

        humann3 --version  > versions.txt
        """
}



//This process merges the taxonomy tables generated by metaphlan
process COMBINE_READ_BASED_PROCESSING_TAXONOMY {

    tag "Merging metaphlan taxonomy tables..."
    label "read_based"
    label "read_based_outputs"

    input:
        path(metaphlan_bugs_list_files)
    output:
        path("${params.additional_filename_prefix}Metaphlan-taxonomy${params.assay_suffix}.tsv"), emit: taxonomy
        path("versions.txt"), emit: version
    script:
        """
        merge_metaphlan_tables.py ${metaphlan_bugs_list_files} \\
                         > ${params.additional_filename_prefix}Metaphlan-taxonomy${params.assay_suffix}.tsv 2> /dev/null

        # Removing redundant text from headers 
        sed -i 's/_metaphlan_bugs_list//g' ${params.additional_filename_prefix}Metaphlan-taxonomy${params.assay_suffix}.tsv

        metaphlan --version > versions.txt
        """
}



workflow read_based {

    take:
        filtered_reads
        krakendb_dir
        kaijudb_dir
        chocophlan_dir
        uniref_dir
        metaphlan_dir
        utilities_dir

    main:
        
        software_versions_ch = Channel.empty()
        // Kraken
        if(krakendb_dir){
            KRAKEN_CLASSIFY(krakendb_dir, filtered_reads)
        }else{
            SETUP_KRAKEN(params.krakendb_url)
            KRAKEN_CLASSIFY(SETUP_KRAKEN.out.krakendb_dir, filtered_reads)
            SETUP_KRAKEN.out.version | mix(software_versions_ch) | set{software_versions_ch}
        }
        kraken_reports = KRAKEN_CLASSIFY.out.report.map{sample_id, report -> report}.collect()
        KRAKEN2TABLE(kraken_reports)
        KRAKEN2KRONA(KRAKEN_CLASSIFY.out.report)
        KRAKEN_REPORT("kraken", KRAKEN2KRONA.out.krona.collect())

        // Kaiju
        if(kaijudb_dir){
           KAIJU_CLASSIFY(kaijudb_dir, filtered_reads)
        }else{
            SETUP_KAIJU(params.kaijudb_name)
            KAIJU_CLASSIFY(SETUP_KAIJU.out.kaijudb_dir, filtered_reads)
            SETUP_KAIJU.out.version | mix(software_versions_ch) | set{software_versions_ch}
        }
        kaiju_reports = KAIJU_CLASSIFY.out.report.map{sample_id, report -> report}.collect()
        KAIJU2TABLE(params.kaijudb_dir, "species", kaiju_reports)
        KAIJU2KRONA(params.kaijudb_dir, KAIJU_CLASSIFY.out.report)
        KAIJU_REPORT("kaiju", KAIJU2KRONA.out.krona.collect())

        // Humann
        if(chocophlan_dir && uniref_dir && metaphlan_dir && utilities_dir){
            HUMANN(filtered_reads, chocophlan_dir, uniref_dir, metaphlan_dir)
        }else{
            make_humann_db()
            HUMANN(filtered_reads, make_humann_db.out.chocophlan_dir, 
                   make_humann_db.out.uniref_dir,
                   make_humann_db.out.metaphlan_db_dir) 
            make_humann_db.out.versions | mix(software_versions_ch) | set{software_versions_ch} 
        }
        METAPHLAN2KRONA(HUMANN.out.metaphlan_bugs_list)
        METAPHLAN_REPORT("metaphlan", METAPHLAN2KRONA.out.krona.collect())

        gene_families_ch = HUMANN.out.genefamilies.collect()
        pathabundance_ch = HUMANN.out.pathabundance.collect()
        pathcoverage_ch = HUMANN.out.pathcoverage.collect()
        metaphlan_bugs_list_ch = HUMANN.out.metaphlan_bugs_list.map{sample_id, bug_list -> bug_list}.collect()
        
        if(chocophlan_dir && uniref_dir && metaphlan_dir && utilities_dir){
            COMBINE_READ_BASED_PROCESSING_TABLES(gene_families_ch, pathabundance_ch, pathcoverage_ch, utilities_dir)
        }else{
            COMBINE_READ_BASED_PROCESSING_TABLES(gene_families_ch, pathabundance_ch,
                                              pathcoverage_ch, make_humann_db.out.utilities_dir)
        }
        SPLIT_READ_BASED_PROCESSING_TABLES(COMBINE_READ_BASED_PROCESSING_TABLES.out.gene_families,
                                           COMBINE_READ_BASED_PROCESSING_TABLES.out.path_abundances,
                                           COMBINE_READ_BASED_PROCESSING_TABLES.out.path_coverages)

        GEN_NORMALIZED_READ_BASED_PROCESSING_TABLES(SPLIT_READ_BASED_PROCESSING_TABLES.out.gene_families,
                                                    SPLIT_READ_BASED_PROCESSING_TABLES.out.path_abundances)

        GEN_READ_BASED_PROCESSING_KO_TABLE(SPLIT_READ_BASED_PROCESSING_TABLES.out.gene_families)
        ko_table_ch = GEN_READ_BASED_PROCESSING_KO_TABLE.out.gene_families
        
        COMBINE_READ_BASED_PROCESSING_TAXONOMY(metaphlan_bugs_list_ch)
        taxonomy_ch = COMBINE_READ_BASED_PROCESSING_TAXONOMY.out.taxonomy

        
        KRAKEN_CLASSIFY.out.version | mix(software_versions_ch) | set{software_versions_ch}
        KRAKEN2TABLE.out.version | mix(software_versions_ch) | set{software_versions_ch}
        KRAKEN2KRONA.out.version | mix(software_versions_ch) | set{software_versions_ch}
        KRAKEN_REPORT.out.version  | mix(software_versions_ch) | set{software_versions_ch}
        KAIJU_CLASSIFY.out.version | mix(software_versions_ch) | set{software_versions_ch}
        KAIJU2TABLE.out.version  | mix(software_versions_ch) | set{software_versions_ch}
        KAIJU2KRONA.out.version  | mix(software_versions_ch) | set{software_versions_ch}
        KAIJU_REPORT.out.version | mix(software_versions_ch) | set{software_versions_ch}
        HUMANN.out.version | mix(software_versions_ch) | set{software_versions_ch}
        METAPHLAN_REPORT.out.version | mix(software_versions_ch) | set{software_versions_ch}
        
        COMBINE_READ_BASED_PROCESSING_TABLES.out.version | mix(software_versions_ch) | set{software_versions_ch}
        SPLIT_READ_BASED_PROCESSING_TABLES.out.version | mix(software_versions_ch) | set{software_versions_ch}
        GEN_NORMALIZED_READ_BASED_PROCESSING_TABLES.out.version | mix(software_versions_ch) | set{software_versions_ch}
        GEN_READ_BASED_PROCESSING_KO_TABLE.out.version | mix(software_versions_ch) | set{software_versions_ch}
        COMBINE_READ_BASED_PROCESSING_TAXONOMY.out.version | mix(software_versions_ch) | set{software_versions_ch}

    emit:
        gene_families = GEN_NORMALIZED_READ_BASED_PROCESSING_TABLES.out.gene_families
        path_abundances = GEN_NORMALIZED_READ_BASED_PROCESSING_TABLES.out.path_abundances
        ko_table = ko_table_ch
        taxonomy = taxonomy_ch 
        versions = software_versions_ch
}


workflow {

    read_based(filtered_ch, 
                params.krakendb_dir,
                params.kaijudb_dir,
                params.chocophlan_dir,
                params.uniref_dir,
                params.metaphlan_db_dir,
                params.utilities_dir)
}
