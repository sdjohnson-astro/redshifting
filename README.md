# redshifting

# Set up environment
export PATH="parentdirectorytree/redshifting:$PATH"
export PYTHONPATH="parentdirectorytree/redshifting/:$PYTHONPATH"
export REDSHIFTING="parentdirectorytree/redshifting"

# Running the GUI:
cd into the directory for a mask reduced with CosmosPipeline
and enter this command into the terminal: cubs_redshifting.py -m maskname


# Bugs:
-- Currently if the user changes the x size of the app the 2D and 1D spectra no longer align
   because the plot margins more relative to one another
   
-- redshifts in table are not being shown to consistent number of digits

# Annoyances
-- zooming recenters on mouse position is annoying me a bit. need input from others