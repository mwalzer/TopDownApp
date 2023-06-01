## Workflows
This folder contains the workflows as used by the app (folder: `local`) as well as different configurations for batch processing (folder: `cluster`). For infrastructure specialisation configuration, please see the `config` folder. The workflows from the `local` folder do not need any extra configuration for container as they are meant to run either locally (with the necessary tools installed) or inside the `tda-app` container (e.g. `tda-app-latest`) and have no container instructions inside the respective processes.

### Examples

#### Run locally like the app
For this, you will need all necessary tools installed, or more easy, start from within the app container:
```
singularity shell containers/topdown-app\:may23.simg 
Singularity> cd /opt/app/
Singularity> nextflow run local/topdown_params.nf --raw_file /home/data/st_1.raw --fasta_file /home/data/tutorial_uniprot-st.fasta --mods modconf/common_mods_Ox.txt 
```

#### Run on HPC
Running any automated workflow on HPC infrastructure will require specialised configuration. For an LSF cluster for example, you need to enable LSF scheduling and singularity, and successively instruct the location of the singularity images as seen from a cluster node, like in `nf_lsf.config` and `cluster_BPA.yml`. A copy of the repository on the cluster simplifies structuring your analysis project and access to necessary (config) files.
```
nextflow run 'workflow/cluster/topdown_cluster_noconversion_FD+PIC.nf' -c 'config/nf_lsf.config' -params-file 'config/cluster_BPA.yml' --mzML_file $FILE --outdir $ODIR
```

### Run with a newly developed workflow
Adding a new tool usually starts with creating a tool container, and then adding a process using this tool to the nextflow at the adequate position. The process can be instructed to use the new container either directly in the nextflow or via a `params-file`.
```
nextflow run workflow/local/topdown_local_MSpathfinderT.nf -params-file config/mspt.yml -c config/nf.config --raw_file /home/data/st_1.raw --fasta_file /home/data/tutorial_uniprot-st.fasta
```

#### Create a mermaid markdown document to overview your workflow 
You can create a visual overview represnetation of your workflow with the topdown.mmd file on GitHub or [mermaid online](https://mermaid.live).
```
nextflow run local/topdown_local.nf -with-dag topdown.mmd --raw_file '/home/data/st_2.raw' --fasta '/home/data/tutorial_uniprot-st.fasta' --mods 'config/common_mods_OxMeth.txt'
``` 