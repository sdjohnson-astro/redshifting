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

-- in python 3 strings are weird in the table. e.g. b'comment'


# Annoyances
-- zooming recenters on mouse position is annoying me a bit. need input from others

-- currently updating table entries by replacing all data in the table. This re-scrolls to the top which is annoying. There must be a better way.

-- percentiles on 2D image are not saved when you change the redshift