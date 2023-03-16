## containers
This folder contains the build recipes for various containers for singularity.
	
build example: `singularity build --fakeroot containers/topdown-main\:dev.simg containers/topdown-main.sdef`

docker(local) to singularity conversion example: 
`singularity build my_container.simg docker-daemon://local/my_container`

docker(hub) to singularity conversion example: 
`singularity build TPP.simg docker://spctools/tpp`