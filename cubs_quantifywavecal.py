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
from scipy.stats import norm

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
def getsky5580(id):
   
   specfilename = '{}{}_1dspec.fits'.format(directory_carpy, id)
   spec = fits.getdata(specfilename)
   header = fits.getheader(specfilename)
   
   wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
   skyflux = spec[2, :]
   
   
   sky = Table()
   sky['wave'] = wave
   sky['sky'] = skyflux
   
   sky = sky[(sky['wave'] > 5550) & (sky['wave'] < 5610)]
   
   cdelt1 = header['CDELT1']
   
   return sky, cdelt1


# Function to get the sky spectrum for a source with input id
def getsky8400(id):
   
   specfilename = '{}{}_1dspec.fits'.format(directory_carpy, id)
   spec = fits.getdata(specfilename)
   header = fits.getheader(specfilename)
   
   wave = header['CRVAL1'] + header['CDELT1']*np.arange(header['NAXIS1'])
   skyflux = spec[2, :]
   
   
   sky = Table()
   sky['wave'] = wave
   sky['sky'] = skyflux
   
   sky = sky[(sky['wave'] > 8377) & (sky['wave'] < 8418)]
   
   cdelt1 = header['CDELT1']
   
   return sky, cdelt1



def fitsky5580(id):
   
   sky, cdelt1 = getsky5580(id)
   
   
   gmodel = Model(gaussian)
   
   parameters = Parameters()
   parameters.add_many(('amp',     np.max(sky['sky']) - np.min(sky['sky']), True, None,    None,    None),
                       ('cen',     5578.5,                   True, 5570.0,  5590.0,  None),
                       ('wid',     5.0,                      True, 1.0,     None,    None),
                       ('a',       0.0,                      True, None,    None,    None),
                       ('b',       np.min(sky['sky']), True, None,    None,    None))
   
   result = gmodel.fit(sky['sky'], params=parameters, x=sky['wave'])
   
   return result, cdelt1

def fitsky8400(id):
   
   sky, cdelt1 = getsky8400(id)
   
   
   gmodel = Model(gaussian)
   
   parameters = Parameters()
   parameters.add_many(('amp',     np.max(sky['sky']) - np.min(sky['sky']), True, None,    None,    None),
                       ('cen',     8401.484,                   True, 8390,  8410,  None),
                       ('wid',     5.0,                      True, 1.0,     None,    None),
                       ('a',       0.0,                      True, None,    None,    None),
                       ('b',       np.min(sky['sky']), True, None,    None,    None))
   
   result = gmodel.fit(sky['sky'], params=parameters, x=sky['wave'])
   print(result.fit_report())
   
   return result, cdelt1

 
def gaussian(x, amp, cen, wid, a, b):
    """1-d gaussian: gaussian(x, amp, cen, wid)"""
    return (amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2)) + a*(x - np.median(x)) + b
 

#objects = objects[0:50]
objects['wave5580'] = np.nan
objects['wave8400'] = np.nan
for object in objects:
   
   # Try the 5580 sky line
   try:
      result, cdelt1 = fitsky5580(object['id'])
      object['wave5580'] = result.best_values['cen']
      
   except:
      print('Object files not found.')
      
   
   # Try the 8400 sky line
   try:
      result, cdelt1 = fitsky8400(object['id'])
      object['wave8400'] = result.best_values['cen']
      
   except:
      print('Object files not found.')

objects = objects[np.isfinite(objects['wave5580']) & np.isfinite(objects['wave8400'])]

print(objects['wave5580'])

print('median = {:0.2f}'.format(np.median(objects['wave5580'])))
print('mean = {:0.2f}'.format(np.mean(objects['wave5580'])))
print('sigma = {:0.2f}'.format(np.std(objects['wave5580'])))
print('resistant sigma = {:0.2f}'.format(median_absolute_deviation(objects['wave5580'])*1.486))
print('dw(obs - vac) = {:0.2f}'.format(np.median(objects['wave5580']) - 5578.5))
print('dw(obs - air) = {:0.2f}'.format(np.median(objects['wave5580']) - 5576.95))


mu_median = np.median(objects['wave5580'])
robust_sigma = median_absolute_deviation(objects['wave5580'])*1.486

fig, ax = plt.subplots(2, figsize=(9, 7))

bins = np.arange(5578.5 - 300.0/c_kms*5578.5, 5578.5 + 600.0/c_kms*5578.5, 0.1)
x = np.arange(5578.5 - 300.0/c_kms*5578.5, 5578.5 + 600.0/c_kms*5578.5, 0.01)
pdf = norm.pdf(x, mu_median, robust_sigma)


ax[0].hist(objects['wave5580'], histtype='step', color='black', label=r'$\rm centroids\ for\ each\ slitlet$', bins=bins, density=True)
ax[0].plot(x, pdf, color='grey', linestyle='--', label=r'$\rm outlier\ resistant\ fit\ median={:0.2f}\ \AA\ \sigma = {:0.2f}\ \AA$'.format(mu_median, robust_sigma))
ax[0].axvline(5578.5, color='blue', label=r'$\rm 5578.5\ (vacuum)\ obs-vac = {:0.2f}\ \AA$'.format(mu_median - 5578.5))
ax[0].axvline(5576.95, color='red', linestyle=':', label=r'$\rm 5576.95\ (air)\ obs-air = {:0.2f}\ \AA$'.format(mu_median - 5576.95))
ax[0].plot([5578.5 - cdelt1/2, 5578.5 + cdelt1/2], [0.5, 0.5], color='blue', linestyle='-.', label=r'$\rm 1\ pixel = {:0.2f}\ \AA$'.format(cdelt1))

ax[0].legend()
ax[0].set_ylabel(r'$\rm frequency$')
ax[0].minorticks_on()



bins = np.arange(8401.48 - 300.0/c_kms*8401.48, 8401.48 + 600.0/c_kms*8401.48, 0.1)
x = np.arange(8401.48 - 300.0/c_kms*8401.48,  8401.48 + 600.0/c_kms*8401.48, 0.01)
pdf = norm.pdf(x, mu_median, robust_sigma)

mu_median = np.median(objects['wave8400'])
robust_sigma = median_absolute_deviation(objects['wave8400'])*1.486
pdf = norm.pdf(x, mu_median, robust_sigma)


ax[1].hist(objects['wave8400'], histtype='step', color='black', label=r'$\rm centroids\ for\ each\ slitlet$', bins=bins, density=True)
ax[1].plot(x, pdf, color='grey', linestyle='--', label=r'$\rm outlier\ resistant\ fit\ median={:0.2f}\ \AA\ \sigma = {:0.2f}\ \AA$'.format(mu_median, robust_sigma))

ax[1].axvline(8401.48, color='blue', label=r'$\rm 8401.48\ (vacuum)\ obs-vac = {:0.2f}\ \AA$'.format(mu_median - 8401.48))
ax[1].axvline(8399.17578, color='red', linestyle=':', label=r'$\rm 8399.18\ (air)\ obs-air = {:0.2f}\ \AA$'.format(mu_median - 8399.18))
ax[1].plot([8401.48 - cdelt1/2, 8401.48 + cdelt1/2], [0.5, 0.5], color='blue', linestyle='-.', label=r'$\rm 1\ pixel = {:0.2f}\ \AA$'.format(cdelt1))
ax[1].minorticks_on()
ax[1].legend()

ax[1].set_ylabel(r'$\rm frequency$')
ax[1].set_xlabel(r'$\rm observed\ wavelength\ [\AA]$')

fig.tight_layout()
plt.savefig('{}_wavecal_scatter.pdf'.format(args.m))


os.system('open {}_wavecal_scatter.pdf'.format(args.m))