#!/usr/bin/env python
import numpy as np
from astropy.table import Table
import argparse

# Set up the command line argument parser
parser = argparse.ArgumentParser(description='Create a region file from a MUSE initial object catalog')
parser.add_argument('-f', metavar='catalog filename', type=str, help='name of the catalog to turned into a region file', required=True)
parser.add_argument('-r', metavar='radius', type=float, help='Radius of the circle to draw', required=False, default=0.6)

args = parser.parse_args()


objects = Table.read('{}.dat'.format(args.f), format='ascii.fixed_width')
print(objects)

reg = open('{}.reg'.format(args.f), 'w')

for object in objects:
   
   reg.write('fk5; circle({}, {}, {}") # text={}\n'.format(object['ra'], object['dec'], args.r, '{' + str(object['id']) + '}'))

reg.close()