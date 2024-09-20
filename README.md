# redshifting



# Dependencies
* numpy
* scipy
* astropy
* PyQt5 (https://www.riverbankcomputing.com/software/pyqt)
* pyqtgraph (http://www.pyqtgraph.org/)
* lmfit (https://lmfit.github.io/lmfit-py/)

# Conda setup
conda create --name redshifting python=3.10.14

conda activate redshifting

conda install -c conda-forge pyqt
conda install -c anaconda pyqtgraph
conda install -c conda-forge astropy
conda install -c conda-forge lmfit
conda install -c astropy corner
conda activate redshifting

# Download the code and add these to your startup file (.zshrc if using conda on a mac)
export REDSHIFTING=“parentdirectorytree/redshifting”
export PATH=“$REDSHIFTING:$PATH”
export PYTHONPATH=“$REDSHIFTING:$PYTHONPATH”
# Running the GUI:
cd into the directory *above* the folder with the reductions for that mask
 (i.e., if you download the ldss3_v3 spectra, the folder you untar them into)
and enter this command into the terminal: cubs_redshifting.py -m maskname

See cubs_redshifting.txt for a usage guide and keystroke listing

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

