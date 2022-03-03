#!/usr/bin/env python
import sys
import numpy as np
from astropy.table import Table, hstack, vstack, join
from matplotlib import pyplot as p
import argparse

c_kms = 299792.458


parser = argparse.ArgumentParser(description='Compared _objects.fits file with consensus table')
parser.add_argument('-m', metavar='mask name', type=str, help='mask name', required=True)
args = parser.parse_args()

objects_filename = '{}_spec1D/{}_objects.fits'.format(args.m, args.m)
objects = Table.read(objects_filename)

print(objects)


concensus = Table.read('{}_concensus.dat'.format(args.m), format='ascii')
print(concensus)


objects = join(objects, concensus, keys=('row'))
objects['dv'] = (objects['redshift'] - objects['redshift_final'])/(1 + objects['redshift_final'])*c_kms
print(objects)

print('************************************************************************')
print('Quality mis-match')
objects_quality = objects[objects['quality'] != objects['quality_final']]
print(objects_quality)


print('************************************************************************')
print('class mis-match')
objects_class = objects[objects['class'] != objects['class_final']]
print(objects_class)


print('************************************************************************')
print('redshift mis-match')
objects_dv = objects[np.abs((objects['dv']) > 20.0) & (objects['quality'] >= 1)]
print(objects_dv)


if (len(objects_quality) == 0) & (len(objects_class) == 0) & (len(objects_dv) == 0):
   
   print('No mismatches found')
   
else:
   print('Mis-matches found!!!! Check please')