# redshifting

# Set up environment

```
export REDSHIFTING="parentdirectorytree/redshifting"
export PATH="$REDSHIFTING:$PATH"
export PYTHONPATH="$REDSHIFTING:$PYTHONPATH"
```

# Dependencies
* numpy
* scipy
* astropy
* PyQt5 (https://www.riverbankcomputing.com/software/pyqt)
* pyqtgraph (http://www.pyqtgraph.org/)
* lmfit (https://lmfit.github.io/lmfit-py/)

# Conda setup

```
conda create -n redshifting
conda install -n redshifting python=3 numpy scipy matplotlib pyqtgraph
conda install -n redshifting astropy
conda install -n redshifting -c GSECARS lmfit
conda install -n redshifting -c anaconda pyqt

source activate redshifting
```


# Running the GUI:
cd into the directory for a mask reduced with CosmosPipeline
and enter this command into the terminal: cubs_redshifting.py -m maskname


# Bugs:
* redshifts in table are not being shown to consistent number of digits
* in python 3 strings are weird in the table. e.g. b'comment'
* histogram for 2D image random triangles show up
* (Old bug with resizing should be fixed)

# Annoyances
* zooming recenters on mouse position is annoying me a bit. need input from others
* currently updating table entries by replacing all data in the table. This re-scrolls to the top which is annoying. There must be a better way. (I fixed this, but I think in a dumb way).
* percentiles on 2D image are not saved when you change the redshift
* implement an advance but no save feature

