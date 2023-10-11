#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// include non-process modules
include { help_message; version_message; complete_message; error_message; pipeline_start_message } from './modules/messages.nf'
include { default_params; check_params } from './modules/params_parser.nf'
include { help_or_version } from './modules/params_utilities.nf'

version = '1.0dev'

// setup default params
default_params = default_params()

// merge defaults with user params
merged_params = default_params + params

// help and version messages
help_or_version(merged_params, version)

final_params = check_params(merged_params)

// starting pipeline
pipeline_start_message(version, final_params)



// include processes
include { MINIMAP2_SAM; SAM_SORT_AND_INDEX; STRAINBERRY; COMBINE_LOG_SUMMARY } from './modules/processes.nf' addParams(final_params)



workflow  {
          read_ch = channel
                          .fromPath( final_params.reads )
                          .map { file -> tuple(file.simpleName, file) }
			  .ifEmpty { error "Cannot find any reads matching: ${final_params.reads}" }

          assemblies_ch = channel
                                .fromPath( final_params.assemblies, checkIfExists: true )
                                .map { file -> tuple(file.simpleName, file) }
				.ifEmpty { error "Cannot find any assemblies matching: ${final_params.assemblies}" }

         joined_ch = read_ch.join(assemblies_ch)

 
	 MINIMAP2_SAM ( joined_ch )
	 
	 joined_sam_ch = MINIMAP2_SAM.out.sam_ch.join(assemblies_ch)

         SAM_SORT_AND_INDEX(joined_sam_ch)
	 
	 joined_sam_assembly_ch = SAM_SORT_AND_INDEX.out.bam_ch.join(assemblies_ch)


         STRAINBERRY(joined_sam_assembly_ch)
	 
	 collected_log_summary_ch = STRAINBERRY.out.log_summary_ch.collect( sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} )

         COMBINE_LOG_SUMMARY(collected_log_summary_ch)
	 

}

workflow.onComplete {
    complete_message(final_params, workflow, version)
}

workflow.onError {
    error_message(workflow)
}
