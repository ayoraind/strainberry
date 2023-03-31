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
    
    
    errorStrategy { task.attempt <= 5 ? "retry" : "finish" }
    maxRetries 5
    
    input:
    tuple val(meta), path(bam), path(bai), path(fai), path(assembly)
    
    output:
    tuple val(meta), path("${meta}.fasta"), emit: first_assembly_ch
    path "versions.yml"                   , emit: versions_ch

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"
    """
    # strain resolution
          
     strainberry $args -r $assembly -b $bam -c $task.cpus -o sberry_out
    
     mv sberry_out/assembly.scaffolds.fa ${prefix}.fasta
    
    sed -i "s/^>/>${prefix}_/g" ${prefix}.fasta
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strainberry: \$( strainberry --version 2>&1 | sed 's/Strainberry //g' )
    END_VERSIONS
    """
}
