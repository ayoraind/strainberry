process MINIMAP2 {
    tag "mapping of $meta assemblies"
    
  //  publishDir "${params.output_dir}", mode:'copy'
    
    
    errorStrategy { task.attempt <= 5 ? "retry" : "finish" }
    maxRetries 5
    
    input:
    tuple val(meta), path(reads), path(assembly)

    output:
    tuple val(meta), path("${meta}.alignment.sorted.bam"),     emit: bam_ch
    path("${meta}.alignment.sorted.bam.bai"),                  emit: bai_ch
    path("${meta}.fasta.fai"),                                 emit: fai_ch
    path "versions.yml",                                       emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    # mapping
     minimap2 -ax map-ont -t $task.cpus $assembly $reads | samtools sort > ${prefix}.alignment.sorted.bam
     samtools faidx $assembly
     samtools index -b ${prefix}.alignment.sorted.bam
     
         
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
         minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}

process MINIMAP2_SAM {
    tag "$meta"

   // publishDir "${params.output_dir}", mode:'copy'

    errorStrategy { task.attempt <= 5 ? "retry" : "finish" }
    maxRetries 5

    input:
    tuple val(meta), path(reads), path(assembly)

    output:
    tuple val(meta), path("*.sam"), emit: sam_ch
    path "versions.yml" , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    minimap2 -ax map-ont -t 12 $assembly $reads  > ${meta}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}

process SAM_SORT_AND_INDEX {
    tag "$meta"

   // publishDir "${params.output_dir}", mode:'copy'

    input:
    tuple val(meta), path(sam), path(assembly)

    output:
    tuple val(meta), path("*.bam"), path("${meta}.alignment.sorted.bam.bai"), path("${meta}.fasta.fai"), emit: bam_ch
    path "versions.yml",                                                                                 emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    """

    samtools sort $sam > ${meta}.alignment.sorted.bam
    samtools faidx $assembly
    samtools index -b ${meta}.alignment.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
	 END_VERSIONS
    """
}

process STRAINBERRY {
    tag "Strain resolution of $meta assemblies"
    
    
    publishDir "${params.output_dir}/${meta}_FLYE_SBERRY", mode:'copy'
    
    
    errorStrategy { task.attempt <= 5 ? "retry" : "ignore" }
    maxRetries 5
    
    input:
    tuple val(meta), path(bam), path(bai), path(fai), path(assembly)
    
    output:
    tuple val(meta), path("${meta}.fasta"), emit: first_assembly_ch
    tuple val(meta), path("${meta}.log"),   emit: log_ch 
    path("*.logsummary"),                   emit: log_summary_ch
    path "versions.yml"                   , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    # strain resolution
          
     strainberry $args -r $assembly -b $bam -c $task.cpus --nanopore -o sberry_out &> ${meta}.log
    
     mv sberry_out/assembly.scaffolds.fa ${prefix}.fasta
    
    sed -i "s/^>/>${prefix}_/g" ${prefix}.fasta
    
    # summarize log file
    log_summary.sh ${meta}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strainberry: \$( strainberry --version 2>&1 | sed 's/Strainberry //g' )
    END_VERSIONS
    """
}


process COMBINE_LOG_SUMMARY {
    publishDir "${params.output_dir}", mode:'copy'
    tag { 'combine log summary files'} 
    
    
    input:
    path(log_summary_files)
    

    output:
    path("combined_log_summary.txt"), emit: combined_log_summary_ch

    
    script:
    """ 
    LOG_SUMMARY_FILES=(${log_summary_files})
    
    for index in \${!LOG_SUMMARY_FILES[@]}; do
    LOG_SUMMARY_FILE=\${LOG_SUMMARY_FILES[\$index]}
    
    # add header line if first file
    if [[ \$index -eq 0 ]]; then
      echo "\$(head -1 \${LOG_SUMMARY_FILE})" >> combined_log_summary.txt
    fi
    echo "\$(awk 'FNR==2 {print}' \${LOG_SUMMARY_FILE})" >> combined_log_summary.txt
    done

    """
}
