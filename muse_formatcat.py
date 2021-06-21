#!/usr/bin/env python
import numpy as np
from astropy.table import Table
import argparse
import coord

# Set up the command line argument parser
parser = argparse.ArgumentParser(description='Create a muse initial input catalog')
parser.add_argument('-f', metavar='catalog filename', type=str, help='name of the catalog to be formatted', required=True)
parser.add_argument('-r', metavar='radius', type=float, help='Initial extraction radius', required=False, default=0.0)
parser.add_argument('-c', metavar='caps', type=bool, help='Initial catalog has caps?', required=False, default=False)

args = parser.parse_args()

objects = Table.read(args.f, format='ascii')

if args.c:
   objects.rename_column('ID', 'id')
   objects.rename_column('RA_J2000', 'ra')
   objects.rename_column('DEC_J2000', 'dec')

print(objects)
ra_string, dec_string, objects['name'] = coord.coordstring(objects['ra'], objects['dec'])

objects['row'] = np.arange(len(objects)) + 1
if args.r == 0:
   objects['radius'] = objects['r_Kron']*0.2
else:
   objects['radius'] = args.r
objects = objects['row', 'id', 'name', 'ra', 'dec', 'radius']


print(objects)

objects.write(args.f, format='ascii.fixed_width', overwrite=True)