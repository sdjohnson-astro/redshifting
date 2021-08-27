#!/usr/bin/env python
import numpy as np
from astropy.table import table
from mpdaf.obj import Cube
import argparse


parser = argparse.ArgumentParser(description='Create a muse whitelight image')
parser.add_argument('-f', metavar='MUSE datacube filename', type=str, help='name of the MUSE to catalog', required=True)
args = parser.parse_args()



cube = Cube(args.f)
whitelight = cube.median(axis=0)
print(whitelight.shape)

savename = args.f
savename = savename.replace('.fits', '_WHITE.fits')
whitelight.write(savename)


#cube = Cube('J0333-4102_COMBINED_CUBE_MED_FINAL_vac_subtracted.fits')
#whitelight = cube.median(axis=0)
#print(whitelight.shape)

#whitelight.write('J0333-4102_COMBINED_SUBTRACTED_WHITE.fits')