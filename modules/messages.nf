def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf --reads "PathToReadFile(s)" --output_dir "PathToOutputDir" --assemblies "PathToAssemblies" 

        Mandatory arguments:
         --reads                        Query fastqz file of sequences you wish to supply as input (e.g., "/MIGE/01_DATA/01_FASTQ/T055-8-*.fastq.gz")
         --assemblies                   Either base assembly (in this case from Flye) or error-corrected assembly should be sufficient (e.g., "/MIGE/01_DATA/03_ASSEMBLY/*_FLYE/T055-8-*.fasta").
         --output_dir	                Output directory

        Optional arguments:
         --help                         This usage statement.
         --version                      Version statement
        """
}

version = '1.0dev'

def Version() {
      println(
            """
            ====================================================
             STRAIN RESOLUTION: TAPIR Pipeline version ${version}
            ====================================================
            """.stripIndent()
        )

}


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Show version
if (params.version) {
    Version()
    exit 0
}


