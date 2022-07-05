#!/usr/bin/env python
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
import argparse
import redshift
import os


# Set up the command line argument parser
parser = argparse.ArgumentParser(description='Perform MCMC likelihood exploration of redshift for a CUBS spectrum')
parser.add_argument('-m', metavar='mask name', type=str, help='mask name', required=True)
parser.add_argument('-r', metavar='row number', type=int, help='row number on mask', required=True, default=False)
parser.add_argument('-z', metavar='Initial guess redshift', type=float, help='Initial guess redshift', required=False, default=-1.0)
parser.add_argument('-nsteps', metavar='number of steps', type=int, help='convert air to vacuum?', required=False, default=2500)
parser.add_argument('-nburn', metavar='burn-in', type=int, help='force convert air to vacuum?', required=False, default=1000)



args = parser.parse_args()



objects = Table.read('{}_spec1D/{}_objects.fits'.format(args.m, args.m))

# find the object and get the spectrum.
object = objects[objects['row'] == args.r]
filename_spec = '{}_spec1D/{}_{}_{}.fits'.format(args.m, args.m, object['row'][0], object['id'][0])
print(filename_spec)
spec = Table.read(filename_spec)

spec = spec[spec['mask'] == 1]

if args.z == -1.0:
   redshift_start = object['redshift'][0]
else:
   redshift_start = args.z


print(object)

print(spec)

print('Redshift for fitting start = {:0.4f}'.format(redshift_start))


filename_corner = filename_spec.replace('.fits', '_corner.pdf')
filename_specplot = filename_spec.replace('.fits', '_spec.pdf')
result_emcee = redshift.fitz_galaxy_emcee(spec, redshift_start, steps=args.nsteps, burn=args.nburn, saveCorner=filename_corner, saveSpecPlot=filename_specplot)

os.system('open {}'.format(filename_corner))
os.system('open {}'.format(filename_specplot))