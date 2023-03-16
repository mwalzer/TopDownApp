#!/usr/bin/env nextflow
nextflow.enable.dsl = 1
/*
=====
            TopDown analysis workflow for PRIDE
=====

 @# Authors
 Mathias Walzer <walzer@ebi.ac.uk>
-----
Pipeline overview:
1. Mzml conversion
2. Apply deconv
3. Use search engine 
4. outputs - original tsvs and a mzTab 
-----
To-Dos:
(nice to have: flashIDA for data producers - for users that do their own data acquisition)
*/

/*
 * Setup of config and vars
 */
params.raw_file = params.raw_file ?: { log.error "No raw file provided. Make sure you have used the '--raw_file' option."; exit 1 }()
params.fasta = params.fasta ?: { log.error "No fasta file provided. Make sure you have used the '--fasta' option."; exit 1 }()
params.mods = params.mods ?: { log.error "No mods file provided. Make sure you have used the '--fasta' option."; exit 1 }()

//params.raw_file = params.raw_file ?: { log.error "No raw file provided. Make sure you have used the '--raw_file' option."; return "/home/walzer/ms-tools/TopDown/TopPIC_tutorial/st_2.raw" }()
//params.fasta = "/home/walzer/ms-tools/TopDown/TopPIC_tutorial/TopPIC_tutorial_uniprot-st.fasta"
//params.mods = "/home/walzer/ms-tools/TopDown/common_mods.txt"
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into '/tmp/results'"; return "/tmp/results" }()
raw_file = file(params.raw_file)

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
    //container '/home/walzer/other-tools/singularity_images/biocontainers-thermorawfileparser-1.2.3.simg'
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

process FLASHDeconv{
    //container '/home/walzer/ms-tools/TopDown/topdown-tools:feb23.simg'
    memory { 4.GB * task.attempt }
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
 
    // errorStrategy { assert task.errorMessage ==~ 'Centroided data provided but profile spectra expected' ? 'ignore' : 'retry' }  //errorReport???
    // errorStrategy 'ignore'
    errorStrategy 'retry'
 
    input:
    file mzML_file from spectra_files_channel.flatten()
 
    // this will be multiple of each type per raw?
    output:
    file "*.msalign" into deconv_spectra_channel_msalign
    file "*.feature" into deconv_spectra_channel_feats
    file "*.mzML" into deconv_spectra_channel_mzmls
	
    // https://www.nextflow.io/docs/latest/channel.html#fromfilepairs
 
    script:
    """
    FLASHDeconv \
       -in $mzML_file \
       -out ${mzML_file.baseName}_fd.tsv \
       -out_topFD ${mzML_file.baseName}.ms1_msalign ${mzML_file.baseName}.msalign \
       -out_topFD_feature ${mzML_file.baseName}.ms1_feature ${mzML_file.baseName}.feature \
       -out_mzml ${mzML_file.baseName}_fd_deconv.mzML \
       -out_annotated_mzml ${mzML_file.baseName}_fd_annot.mzML
    """
}

/*
 * Generate identifications from the deconvolved MS2 spectra+feature
 * Note: due to TopPICs implicit read of {msalign input file basename}.feature 
 * 	         this node needs to do seemingly superfluous steps
 */
process TopPIC{
    //container '/home/walzer/ms-tools/TopDown/topdown-tools:feb23.simg'
    memory { 4.GB * task.attempt }
    
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    // errorStrategy { assert task.errorMessage ==~ 'Centroided data provided but profile spectra expected' ? 'ignore' : 'retry' }  //errorReport???
    // errorStrategy 'ignore'
    errorStrategy 'retry'
 
    input:
    file msalign_file from deconv_spectra_channel_msalign.flatten()
    file msfeat_file from deconv_spectra_channel_feats.flatten()
	// the tiniest issue with the feature file will make toppic dump core (say, like the ms1 features are contained in there too ...)

    output:
    file "*_proteoform_single.tsv" into id_prtf_channel
    file "*_prsm_single.tsv" into id_prsm_channel
 
    script:
    """
    toppic -d -t FDR -T FDR -u 2 -g -i ${params.mods} ${params.fasta} $msalign_file
    """
	// unless both file variables are mentioned and used the files don't get staged 
}

/*
 * Export TopPIC results to mzTab
 */
process mzTab{
    //container '/home/walzer/ms-tools/TopDown/topdown-tools:feb23.simg'
    memory { 4.GB * task.attempt }
    
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    // errorStrategy { assert task.errorMessage ==~ 'Centroided data provided but profile spectra expected' ? 'ignore' : 'retry' }  //errorReport???
    // errorStrategy 'ignore'
    errorStrategy 'retry'
 
    input:
    file proteoform_ids from id_prtf_channel.flatten()
    file prsm_ids from id_prsm_channel.flatten()
    
    output:
	file "*.mzTab" into mztab_channel
    file "*.h5" into df_channel

    //needs v4 of export_mztab
    script:
    """
    python3 /opt/app/export_mztab.py -s ${params.fasta} -p ${prsm_ids} -f ${proteoform_ids} -5 ${prsm_ids.baseName}.h5 ${prsm_ids.baseName}.mzTab
	"""
}




