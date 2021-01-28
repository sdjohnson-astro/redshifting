#!/usr/bin/env python
import glob
import argparse
from astropy.table import Table
import numpy as np

# Set up the command line argument parser
parser = argparse.ArgumentParser(description='Compare two versions of spec1D files from CUBS IMACS or LDSS3')
parser.add_argument('-d1', metavar='directory 1', type=str, help='Parent directory 1', required=True)
parser.add_argument('-d2', metavar='directory 2', type=str, help='Parent directory 2', required=True)
parser.add_argument('-m', metavar='maskname', type=str, help='mask name', required=True)

args = parser.parse_args()


mask = Table.read('{}/{}_spec1D/{}_objects.fits'.format(args.d1, args.m, args.m))
mask['maxabsDflux'] = 0.0

for object in mask:
   
   try:
   
      filename1 = '{}/{}_spec1D/{}_{}_{}.fits'.format(args.d1, args.m, args.m, object['row'], object['id'])
      spec1 = Table.read(filename1)
      
      filename2 = '{}/{}_spec1D/{}_{}_{}.fits'.format(args.d2, args.m, args.m, object['row'], object['id'])
      spec2 = Table.read(filename2)
      
      print(np.max(np.abs(spec1['flux'] - spec2['flux'])))
      object['maxabsDflux'] = np.max(np.abs(spec1['flux'] - spec2['flux']))
   
   except:
      
      print('file not found')
   
print(mask)

maxabsDiff = np.max(mask['maxabsDflux'])

if maxabsDiff > 0.0:
   
   print('Differences found!!!!!!!!!!!')
   
else:
   
   print('No difference -- ok')