process MINIMAP2 {
    tag "mapping of $meta assemblies"
    
    
    publishDir "${params.output_dir}", mode:'copy'
    
    conda '../sberry_env.yml'
    
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


process STRAINBERRY {
    tag "Strain resolution of $meta assemblies"
    
    
    publishDir "${params.output_dir}/${meta}_FLYE_SBERRY", mode:'copy'
    
    conda '../sberry_env.yml'
    
    input:
    tuple val(meta), path(assembly), path(bam)
    path bai
    path fai
    
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
