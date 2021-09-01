#!/usr/bin/env python
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
import argparse
import glob
from lmfit import Model, Parameters
import os
from pydl.goddard.astro import airtovac
from astropy import units as u
from astropy.io import fits
from astropy.stats import median_absolute_deviation

dv_vacair = 82.68 # km/s at 7000 Ang. At 4000 and 10000 it is 84.76 and 82.46 km/s respectively so constant dv is a very good approximation.
c_kms = 299792.458


# Set up the command line argument parser
parser = argparse.ArgumentParser(description='Correct wavelength and redshifts for air-vacuum issue.')
parser.add_argument('-m', metavar='mask name', type=str, help='name of the mask to be summarized', required=True)


args = parser.parse_args()


directory_spec1D = './{}_spec1D/'.format(args.m)
directory_carpy = './{}/'.format(args.m)

objects = Table.read(directory_spec1D + '{}_objects.fits'.format(args.m))

objects = objects[objects['alignbox'] == False]
print(objects)



# Function to get the sky spectrum for a source with input id
def getsky(id):
   
   specfilename = '{}{}_1dspec.fits'.format(directory_carpy, id)
   spec = fits.getdata(specfilename)
   header = fits.getheader(specfilename)
   
   wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
   skyflux = spec[2, :]
   
   
   sky = Table()
   sky['wave'] = wave
   sky['sky'] = skyflux
   
   sky = sky[(sky['wave'] > 5550) & (sky['wave'] < 5610)]
   
   return sky
 
def gaussian(x, amp, cen, wid, a, b):
    """1-d gaussian: gaussian(x, amp, cen, wid)"""
    return (amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2)) + a*x + b
 

#objects = objects[0:50]
objects['wave5580'] = np.nan
for object in objects:
   
   try:
      sky = getsky(object['id'])
      
      
      gmodel = Model(gaussian)
      
      parameters = Parameters()
      parameters.add_many(('amp',     np.max(sky['sky']) - np.min(sky['sky']), True, None,    None,    None),
                          ('cen',     5578.5,                   True, 5570.0,  5590.0,  None),
                          ('wid',     5.0,                      True, 1.0,     None,    None),
                          ('a',       0.0,                      True, None,    None,    None),
                          ('b',       np.min(sky['sky']), True, None,    None,    None))
      
      result = gmodel.fit(sky['sky'], params=parameters, x=sky['wave'])
      object['wave5580'] = result.best_values['cen']
      print(result.fit_report())
   except:
      print('Object files not found.')

objects = objects[np.isfinite(objects['wave5580'])]

print(objects['wave5580'])

print('median = {:0.2f}'.format(np.median(objects['wave5580'])))
print('mean = {:0.2f}'.format(np.mean(objects['wave5580'])))
print('sigma = {:0.2f}'.format(np.std(objects['wave5580'])))
print('resistant sigma = {:0.2f}'.format(median_absolute_deviation(objects['wave5580'])*1.486))
print('dw(obs - vac) = {:0.2f}'.format(np.median(objects['wave5580']) - 5578.5))
print('dw(obs - air) = {:0.2f}'.format(np.median(objects['wave5580']) - 5576.95))