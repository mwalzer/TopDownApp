# Quick start
Each step involves the use of the <ins>command line</ins>, for _windows_ please start __PowerShell__, on _linux and mac_ open a __terminal__. 

1. Install container engine (or check if it is installed):  

    - Docker: There are binary installation packages available for all platforms: https://docs.docker.com/engine/install/ 
    - Podman: 'drop-in' replacement for Docker: https://podman-desktop.io/downloads 
    - Singularity: https://apptainer.org/admin-docs/master/installation.html 

    To test the installation, type `docker`/`podman`/`singularity` in the command line. If the container engine is installed, you should receive software version and basic operation description printed on screen and no errors. 

2. "Pull" the image: 

    - Docker: `docker pull quay.io/mwalzer/topdownapp:tda-app-latest` 
    - Podman: `podman pull quay.io/mwalzer/topdownapp:tda-app-latest` 
    - Singularity: `singularity build --fakeroot topdownapp.simg docker://quay.io/mwalzer/topdownapp:tda-app-latest` 
    
    (Note that the singularity command creates a a `.simg` image file, so keep in mind to specify the path to the file when "running" the container. However, in singularity, you can "run" docker images on the fly without the image file) 

3. "Run" the image: 
    Since you will operate TopDownApp through the browser, you will have to specify a port mapping –p <from>:<to>, where the from part is defined in TopDownApp by default as port 5006. Changing this will influence the browser address you need to navigate to.   

    - Docker: `docker run -p 5006:5006 quay.io/mwalzer/topdownapp:tda-app-latest` 
    - Podman: `podman run -p 5006:5006 quay.io/mwalzer/topdownapp:tda-app-latest` 
    - Singularity: `singularity run topdownapp.simg` 
    
    (No port mapping necccsesary in singularity. To run docker images in singularity on the fly you can specify a docker registry instead of the image, e.g.: `singularity run docker://quay.io/mwalzer/topdownapp:tda-app-latest`

4. Stop the container 

    - Docker: `CTRL+C` 
    - Podman: `CTRL+C` 
    - Singularity: `CTRL+C` 

5. Mount your data: 
    In certain environments it might be necessary to mount the file system folder which contains your data (<from>) into the container (<to> at `/opt/data`). For example, for `D:\test\` use `-v /d/test:/opt/data`  . 

    - Docker: `docker run -p 5006:5006 -v <from>:<to> quay.io/mwalzer/topdownapp:tda-app-latest` 
    - Podman: `podman run -p 5006:5006 -v <from>:<to> quay.io/mwalzer/topdownapp:tda-app-latest` 
    - Singularity: `singularity run -B <from>:<to> docker://quay.io/mwalzer/topdownapp:tda-app-latest` 

6. Use your browser to navigate to `http://localhost:5006/topdownvisapp`  
