## containers
This folder contains the build recipes for various containers for singularity.
	
build example: `singularity build --fakeroot containers/topdown-main\:dev.simg containers/topdown-main.sdef`

docker(local) to singularity conversion example: 
`singularity build my_container.simg docker-daemon://local/my_container`

docker(hub) to singularity conversion example: 
`singularity build TPP.simg docker://spctools/tpp`


### in stand-alone mode 
For switching parameter setting files, tool versions, etc. it is easier to use the main container that expects a container configuration file to run. helpful also for cluster deployment of the workflow. 

build app ontop tools ontop raw


container depend on each other 
singularity build --fakeroot topdown-tools\:may23.simg topdown-tools.sdef
