#!/usr/bin/env bash
export JAVA_HOME=/hps/software/users/juan/pride/java/jdk-11.0.2/
export PATH=$PATH:$JAVA_HOME
FILE=$1
BASEFILE=$(basename -- "$FILE")
BASE=${FILE%%.*}
mkdir -p /hps/nobackup/juan/pride/topdown/PXD026159_uniprot-sprot_9606_OxMeth
mkdir -p /hps/nobackup/juan/pride/topdown/PXD026159_uniprot-sprot_9606_OxMeth/${BASE}
cd /hps/nobackup/juan/pride/topdown/PXD026159_uniprot-sprot_9606_OxMeth/${BASE}
/hps/nobackup/juan/pride/topdown/nextflow_22.04.5 run '/hps/nobackup/juan/pride/topdown/topdown_cluster_noconversion.nf' -c '/hps/nobackup/juan/pride/topdown/nf.config' --mzML_file $FILE -params-file '/hps/nobackup/juan/pride/topdown/BPA-reident/PXD026159_uniprot-sprot_9606_OxMeth.yml'

