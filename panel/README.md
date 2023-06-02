## Python panel
This folder contains the main app files for panel control, io management, and visualisation.


### Run the app
You can run the app through the container like so:
```
singularity run containers/topdown-main\:dev.simg 
```
Or separately, however you need to make sure the following:
* requirements as listed in `requirements.txt` are met
* all tools used in the workflow process list are installed
* alternatively start from the app container:
```
singularity shell containers/topdown-main\:dev.simg 
Singularity> panel serve panel/topdownvisapp.py
```
And start panel like above. You can also start a local copy of the app while using the app container to cover all dependencies (i.e. not necessarily `/opt/app/topdownvisapp.py`).

For developers the autoreload feature may be of use:
```
panel serve panel/topdownvisapp.py --autoreload
```

### Structure
The app's code is split into three parts:
1. The panel code itself, containing logic and layout - [topdownvisapp](topdownvisapp.py)
2. Data I/O - [ms_io_utils](ms_io_utils.py)
3. Visualisation - [panel_plot_utils](panel_plot_utils.py)


All dependencies are listed in `requirements.txt`. 2. and 3. can be directly reused, e.g. in ipython notebooks, see there. 1. controls the layout of the app, managing data input, and the internal scheduling of the nextflow, and the logic of interactive elments in the anlysis results presentation.

