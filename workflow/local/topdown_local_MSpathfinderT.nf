#!/usr/bin/env nextflow
nextflow.enable.dsl = 1
/*
=====
            TopDown analysis workflows for PRIDE
=====
*/

/*
 * Setup of config and vars
 */
params.raw_file = params.raw_file ?: { log.error "No mzML file provided. Make sure you have used the '--raw_file' option."; exit 1 }()
params.fasta_file = params.fasta_file ?: { log.error "No fasta file provided. Make sure you have used the '--fasta_file' option."; exit 1 }()
//params.mods_file = params.mods_file ?: { log.error "No mods file provided. Make sure you have used the '--mods_file' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()
raw_file = file(params.raw_file)
fasta_file = file(params.fasta_file)

def helpMessage() {
   log.info"""
   ======
   Usage:
   ...

   """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
   helpMessage()
   exit 0
}


/*
 * Generate the mzML + metadata for each RAW file
 */
process ConvertRaw {
    container "${params.allinone.container}"
    memory { 4.GB * task.attempt }
    errorStrategy 'retry'
     
    /*publishDir "${params.outdir}", mode: 'copy', overwrite: true*/
  
    input:
    //file rawFile from "${params.raw_file}"
    
    output: 
    file '*.mzML' into spectra_files_channel
 
    script:
    """
	echo ${params.raw_file}
    ThermoRawFileParser.sh -i=$raw_file -f=2 -o=./
    """
}


/*
 * Generate identifications from the deconvolved MS2 spectra+feature
 */
process MSpathfinderT{
    container "${params.MSpathfinderT.container}"
    memory { 4.GB * task.attempt }

    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    // errorStrategy { assert task.errorMessage ==~ 'Centroided data provided but profile spectra expected' ? 'ignore' : 'retry' }  //errorReport???
    // errorStrategy 'ignore'
    errorStrategy 'retry'

    input:
    file mzML_file from spectra_files_channel.flatten()

    output:
    file "_IcTarget.tsv" into id_prsm_channel
    file "*.mzid" into mzid_channel

    script:
    """
    MSPathFinderT.exe -d .${fasta_file} -s $mzML_file
    """
}
