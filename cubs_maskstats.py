#!/usr/bin/env python
import numpy as np
from astropy.table import Table
import argparse


# Set up the command line argument parser
parser = argparse.ArgumentParser(description='Basic statistics on a mask')
parser.add_argument('-m', metavar='mask name', type=str, help='name of the mask to be summarized', required=True)
args = parser.parse_args()

filename = '{}_spec1D/{}_objects.fits'.format(args.m, args.m)
print('filename = {}'.format(filename))
slits = Table.read(filename)

nSlits = len(slits)

# Good galaxies
galaxies = slits[(slits['class'] == 'galaxy') & (slits['quality'] == 2)]
nGalaxies = len(galaxies)

# Good stars
stars = slits[(slits['class'] == 'star') & (slits['quality'] == 2)]
nStars = len(stars)

# Good quasars
quasars = slits[(slits['class'] == 'quasar') & (slits['quality'] == 2)]
nQuasars = len(quasars)

# Single-line anything
singles = slits[(slits['quality'] == 1)]
nSingles = len(singles)

# failures
failures = slits[(slits['quality'] <= 0)]
nFailures = len(failures)

# succeses rate
successRate = float(nGalaxies + nStars + nQuasars)/float(nSlits)

print('nSlits      = {}'.format(nSlits))
print('nGalaxies   = {}'.format(nGalaxies))
print('nStars      = {}'.format(nStars))
print('nQuasars    = {}'.format(nQuasars))
print('nSuccesses  = {}'.format(nGalaxies + nStars + nQuasars))
print('nSingles    = {}'.format(nSingles))
print('nFailures   = {}'.format(nFailures))
print('successes   = {:0.2f}'.format(successRate))