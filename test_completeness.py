#!/usr/bin/env python
from astropy.io import fits
import argparse
import numpy as np
from astropy.table import Table


parser = argparse.ArgumentParser(description='Verify foreground/background quasars and measure absorbing gas around the foreground in the background quasar spectrum')
parser.add_argument('-m', metavar='mask name', type=str, help='name of the mask to be redshifted', required=True)
args = parser.parse_args()
objects = Table.read('{0}_spec1D/{0}_objects.fits'.format(args.m))
objects=objects[~objects['alignbox']]

ntot=len(objects)
n2=sum(objects['quality']==2)
n2gal=sum((objects['quality']==2) & ((objects['class']=='galaxy') | (objects['class']=='quasar')))
n2star=sum((objects['quality']==2) & (objects['class']=='star'))

n1=sum(objects['quality']==1)
n0=sum(objects['quality']==0)
nf=sum(objects['quality']==-1)

print '{:d} targets\n{:d} Q=2, {:d} galaxies, {:d} stars\n{:d} single line\n{:d} No redshift\n{:d} Bad spectrum'.format(ntot,n2,n2gal,n2star,n1,n0,nf)



