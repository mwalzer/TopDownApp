[![Docker Repository on Quay](https://img.shields.io/badge/container-ready-brightgreen.svg?style=for-the-badge&logo=appveyor "Docker Repository on Quay")](https://quay.io/repository/mwalzer/topdownapp?tab=tags)



# TopDownApp
This repository contains the resources and documentation to develop and build the TopDownApp.
The different components of the app can be found in the respective subfolders.

![Pipeline overview](TopDownAppOverview.png)

The app lets you analyse (1. deconvolve, 2. identify, 3. visualise) TopDown data without having to install multiple software, and is available on every platform.
In the default GUI setting, for the analysis it is expected as input a TopDown LC-MS/MS (for the time being, Thermo) __RAW file__ and matching __FASTA file__ for sequence based identification.
TopDownApp also works without a GUI, e.g. for automated analysis on HPC infrastructure.
It is a highly configurable platform intended to share analysis capability fast, through flexible code reuse. 

## Quick Start
The fastest way to get started with TopDownApp is if you have a container engine already installed ([Singularity](https://apptainer.org/admin-docs/master/installation.html)or [Docker](https://docs.docker.com/engine/install/)).
Then, you can simply start the app right away, then follow the localhost link in the output, usually 

[http://localhost:5006/topdownvisapp](http://localhost:5006/topdownvisapp)!

(See container engine specific commands below.)

### Singularity
It is fastest for Singularity users to pull a pre-made container image from quay.io (docker), due to the lack of a free and open singularity image registry.
```
singularity run -B <e.g. ~/TopDown/data/>:</opt/data/> docker://quay.io/mwalzer/topdownapp:tda-app-latest
```

### Docker
```
docker run -p 5006:5006 -v <e.g. ~/TopDown/data/>:</opt/data/> quay.io/mwalzer/topdownapp:tda-app-latest 
```

For more details on container use, see the `containers` folder.
Or go explore the repository resources for reuse and development below.

## Workflow 
The TopDownApp runs nextflow workflows to effect the analysis of data input. 
These can be reused for data analysis on local compute infrastructure.
See [that folder's README](workflow/README.md).

## Containers
The workflows use containerised tools and the TopDownApp itself can be containerised. 
All build recipes (.sdef) can be found here. 
See [that folder's README](containers/README.md).

## mzTab
The TopDownApp handles input and output data with HUPO-PSI standard formats. 
Development code, dedicated im-/exporters, and documentation can be found here. 
See [that folder's README](mzTab/README.md).

## Python panel
The TopDownApp uses panel to build the browser interface for the app. 
This is also where the core of the TopDownApp lives. 
See [that folder's README](panel/README.md). 

## Notebooks
This folder contains notebooks illustrating TopDownApp re-use/development and (user-submitted) use-case studies.
See [that folder's README](notebooks/README.md). 
