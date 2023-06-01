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
1. Apply deconvolution
2. Use search engine
3. Output - original analysis tsvs and a mzTab
-----
Note that this workflow needs compatible mzML files as input (pwiz conversion for TopSuite)
Note also that this workflow assumes the fasta used to be top indexed already before running 
parallel executions (e.g. by running the workflow once before batch submission)
*/

/*
 * Setup of config and vars
 */
params.mzML_file = params.mzML_file ?: { log.error "No mzML file provided. Make sure you have used the '--mzML_file' option."; exit 1 }()
params.fasta_file = params.fasta_file ?: { log.error "No fasta file provided. Make sure you have used the '--fasta_file' option."; exit 1 }()
params.mods_file = params.mods_file ?: { log.error "No mods file provided. Make sure you have used the '--mods_file' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()
mzML_file = file(params.mzML_file)
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
 * SKIP Generate the mzML + metadata for each RAW file
 */

process FLASHDeconv{
    container "${params.allinone.container}"
    memory { 4.GB * task.attempt }
    publishDir "${params.outdir}", mode: 'copy', overwrite: true
 
    // errorStrategy { assert task.errorMessage ==~ 'Centroided data provided but profile spectra expected' ? 'ignore' : 'retry' }  //errorReport???
    // errorStrategy 'ignore'
    errorStrategy 'retry'
 
    input:
    // file mzML_file from spectra_files_channel.flatten()

    output:
    file "*_ms2.msalign" into deconv_spectra_channel_msalign
    file "*_ms1.feature" into deconv_spectra_channel_feats_ms1
    file "*_ms2.feature" into deconv_spectra_channel_feats_ms2
    file "*.mzML" into deconv_spectra_channel_mzmls
 
    script:
    """
    FLASHDeconv \
       -in $mzML_file \
       -out ${mzML_file.baseName}.tsv \
       -out_topFD ${mzML_file.baseName}_ms1.msalign ${mzML_file.baseName}_ms2.msalign \
       -out_topFD_feature ${mzML_file.baseName}_ms1.feature ${mzML_file.baseName}_ms2.feature \
       -out_mzml ${mzML_file.baseName}_deconv.mzML \
       -out_annotated_mzml ${mzML_file.baseName}_annot.mzML
    """
}

/*
 * Generate identifications from the deconvolved MS2 spectra+feature
 */
process TopMG{
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
    file "*_proteoform.tsv" into allacc_prtf_channel
    file "*_prsm.tsv" into allacc_prsm_channel

    script:
    """
    topmg -d -t FDR -T FDR -u 2 -g -i ${params.mods_file} ${params.fasta_file} $msalign_file
    """
    // unless both file variables are mentioned and used the files don't get staged 
    // probably also something like # echo ${fasta_file.baseName} | cut -c1-25 as TopPic is sensitive to long filenames
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
    python3 /opt/app/export_mztab.py -s ${params.fasta_file} -p ${prsm_ids} -f ${proteoform_ids} -5 ${prsm_ids.baseName}.h5 ${prsm_ids.baseName}.mzTab
    """
}
