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
2. Apply deconvolution
3. Use search engine 
4. Output - original analysis tsvs and a mzTab 
-----
To-Dos:
(nice to have: flashIDA for data producers - for users that do their own data acquisition)
*/

/*
 * Setup of config and vars
 */
params.mzML_file = params.mzML_file ?: { log.error "No mzML file provided. Make sure you have used the '--mzML_file' option."; exit 1 }()
params.fasta = params.fasta ?: { log.error "No fasta file provided. Make sure you have used the '--fasta' option."; exit 1 }()
params.mods = params.mods ?: { log.error "No mods file provided. Make sure you have used the '--mods' option."; exit 1 }()

//params.mzML_file = params.mzML_file ?: { log.error "No raw file provided. Make sure you have used the '--mzML_file' option."; return "/home/walzer/ms-tools/TopDown/TopPIC_tutorial/st_2.raw" }()
//params.fasta = "/home/walzer/ms-tools/TopDown/TopPIC_tutorial/TopPIC_tutorial_uniprot-st.fasta"
//params.mods = "/home/walzer/ms-tools/TopDown/common_mods.txt"
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()
mzML_file = file(params.mzML_file)

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
 * SKIP Generate the mzML + metadata for each RAW file
 */

process TopFD{
    container "${params.allinone.container}"
    memory { 4.GB * task.attempt }
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
 
    // errorStrategy { assert task.errorMessage ==~ 'Centroided data provided but profile spectra expected' ? 'ignore' : 'retry' }  //errorReport???
    // errorStrategy 'ignore'
    errorStrategy 'retry'
 
    input:
    // file mzML_file from spectra_files_channel.flatten()
 
    output:
    file "*.msalign" into deconv_spectra_channel_msalign
    file "*_ms1.feature" into deconv_spectra_channel_feats_ms1
    file "*_ms2.feature" into deconv_spectra_channel_feats_ms2

    // since TopFD tries to write it's results into the same folder as the input,
    // we need some workaround
    script:
    """
    cp -L $params.mzML_file ./${mzML_file.baseName}.mzML
    topfd -g ./${mzML_file.baseName}.mzML
    """
}

/*
 * Generate identifications from the deconvolved MS2 spectra+feature
 * Note: due to TopPICs implicit read of {msalign input file basename}.feature 
 * 	         this node needs to do seemingly superfluous steps
 */
process TopPIC{
    container "${params.allinone.container}"
    memory { 8.GB * task.attempt }
    
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
    // errorStrategy { assert task.errorMessage ==~ 'Centroided data provided but profile spectra expected' ? 'ignore' : 'retry' }  //errorReport???
    // errorStrategy 'ignore'
    errorStrategy 'retry'
 
    input:
    file msalign_file from deconv_spectra_channel_msalign.flatten()
    file ms1feat_file from deconv_spectra_channel_feats_ms1.flatten()
    file ms2feat_file from deconv_spectra_channel_feats_ms2.flatten()
	// the tiniest issue with the feature file will make toppic dump core (say, like the ms1 features are contained in there too ...)

    output:
    file "*_proteoform_single.tsv" into id_prtf_channel
    file "*_prsm_single.tsv" into id_prsm_channel
 
    script:
    """
    #ls $ms1feat_file $ms2feat_file
    toppic -d -t FDR -T FDR -u 2 -g -i ${params.mods} ${params.fasta} $msalign_file
    """
    // unless both file variables are mentioned and used the files don't get staged 
}

/*
 * Export TopPIC results to mzTab
 */
process mzTab{
    container "${params.allinone.container}"
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




