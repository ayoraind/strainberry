## Workflow to resolve strains from long-read assemblies using Strainberry.
### Usage

```

=================================================================================
 STRAIN RESOLUTION OF ASSEMBLIES USING STRAINBERRY: TAPIR Pipeline version 1.0dev
=================================================================================
 The typical command for running the pipeline is as follows:
        nextflow run main.nf --reads "PathToReadFile(s)" --output_dir "PathToOutputDir" --assemblies "PathToAssemblies"

        Mandatory arguments:
         --reads                        Query fastqz file of sequences you wish to supply as input (e.g., "/MIGE/01_DATA/01_FASTQ/T055-8-*.fastq.gz")
         --output_dir                   Output directory (e.g., "/MIGE/01_DATA/03_ASSEMBLY")
	 --assemblies                   must be Flye assembly (e.g., "/MIGE/01_DATA/03_ASSEMBLY/*_FLYE/T055-8-*.fasta")
         
        Optional arguments:
         --help                         This usage statement.
         --version                      Version statement

```


## Introduction
This pipeline attempts to resolve strains from long-read assemblies. This Nextflow pipeline was adapted from the original author's [github page](https://github.com/rvicedomini/strainberry).  


## Sample command
An example of a command to run this pipeline is:

```
nextflow run main.nf --reads "Sample_files/*.fastq.gz" --output_dir "test2" --assemblies "*.fasta"
```

## Word of Note
This is an ongoing project at the Microbial Genome Analysis Group, Institute for Infection Prevention and Hospital Epidemiology, Üniversitätsklinikum, Freiburg. The project is funded by BMBF, Germany, and is led by [Dr. Sandra Reuter](https://www.uniklinik-freiburg.de/iuk-en/infection-prevention-and-hospital-epidemiology/research-group-reuter.html).


## Authors and acknowledgment
The TAPIR (Track Acquisition of Pathogens In Real-time) team.
