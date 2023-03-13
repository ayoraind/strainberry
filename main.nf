#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// include definitions
include  { helpMessage; Version } from './modules/messages.nf'

// include processes
include { MINIMAP2; STRAINBERRY } from './modules/processes.nf'

log.info """\
    STRAIN RESOLUTION  - TAPIR P I P E L I N E
    ==========================================
    output_dir       : ${params.output_dir}
    """
    .stripIndent()


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
