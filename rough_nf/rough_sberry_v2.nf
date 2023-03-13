/*
 * pipeline input parameters
 */
params.assemblies = "$projectDir/data/ggal/gut.fasta"

params.reads = "$projectDir/data/ggal/gut_{1,2}.fastq"
params.output_dir = "$PWD/results"
nextflow.enable.dsl=2


log.info """\
    STRAIN RESOLUTION  - TAPIR   P I P E L I N E
    ============================================
    output_dir       : ${params.output_dir}
    """
    .stripIndent()


process MINIMAP2 {
    tag "mapping of $meta assemblies"
    
    
    publishDir "${params.output_dir}", mode:'copy'
    
    conda './sberry_env.yml'
    
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
    
    conda './sberry_env.yml'
    
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

workflow  {
          read_ch = channel
                          .fromPath( params.reads, checkIfExists: true )
                          .map { file -> tuple(file.simpleName, file) }
			  
	  assemblies_ch = channel
                                .fromPath( params.assemblies, checkIfExists: true )
				.map { file -> tuple(file.simpleName, file) }
				
	 joined_ch = read_ch.join(assemblies_ch)
	 
	          
         MINIMAP2(joined_ch)
	 
	 joined_map_ch = assemblies_ch.join(MINIMAP2.out.bam_ch)
	 
         STRAINBERRY(joined_map_ch, MINIMAP2.out.bai_ch, MINIMAP2.out.fai_ch)

}
