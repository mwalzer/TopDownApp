## Containers
This folder contains the build recipes for various containers for [Singularity](https://apptainer.org/admin-docs/master/installation.html) and [Docker](https://docs.docker.com/engine/install/).
You can find prebuilt Docker containers at [quay.io/mwalzer/topdownapp](https://quay.io/repository/mwalzer/topdownapp)
You will need a running installation of one of the container engines listed in order to run or build the containers as described in the following. For installation guidance please follow the links in the container engine listing above.

### Build Containers

Note: You will need to change to the root directory of the repository to build the containers as shown here.

#### To build for __Singularity__:

* build example: `singularity build --fakeroot containers/topdown-app.simg containers/topdown-app.sdef`
* docker(local) to singularity conversion example: 
`singularity build <chosen image name>.simg docker-daemon://local/mycontainer`
* docker(hub) to singularity conversion example (here for the TPP suite container): 
`singularity build <chosen image name>.simg docker://spctools/tpp`

#### To build for __Docker__:

build example: `docker build -t topdownapp:tda-app-latest -f containers/topdown-app.docker .`

### Run in Stand-alone Mode 
To simply run the TopDownApp, it is easiest to run the latest app container (tag tda-app-latest):
`docker run quay.io/mwalzer/topdownapp:tda-app-latest`

### Run in Development Mode
For switching parameter setting files, tool versions, etc. it is simplest to also use the latest app container (tag tda-app-latest) but start into a command shell:
`docker run -it quay.io/mwalzer/topdownapp:tda-app-latest bash`
This will give you access to all workflows, configuration files (in `/opt/app/) and the Python Panel app, which you can start like so:
`panel serve /opt/app/topdownvisapp.py`

### Structure
The app container (tag `tda-app-latest`, container definition file `topdown-app`) is built ontop of the tools container (tag `tda-tools-latest`, container definition file `topdown-tools`) that is built ontop of the rawfile conversion container (tag `tda-raw-latest`, container definition file `topdown-raw`). This is to make all essential tools and scripts available from one container while keeping individual development contaienrs independent. As such, for cluster deployment, amongst others, this folder contains also individual tools container definitions.
