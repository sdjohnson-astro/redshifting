#!/usr/bin/env python
import numpy as np
from astropy.table import Table
import argparse
import coord

# Set up the command line argument parser
parser = argparse.ArgumentParser(description='Create a muse initial input catalog')
parser.add_argument('-f', metavar='catalog filename', type=str, help='name of the catalog to be formatted', required=True)
args = parser.parse_args()

objects = Table.read(args.f, format='ascii')

ra_string, dec_string, objects['id'] = coord.coordstring(objects['ra'], objects['dec'])

objects['row'] = np.arange(len(objects))

objects = objects['row', 'id', 'ra', 'dec', 'radius']

print(objects)

objects.write(args.f, format='ascii.fixed_width', overwrite=True)