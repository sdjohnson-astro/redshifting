#!/usr/bin/env python
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
import argparse
import glob
from lmfit import Model


def getspecfilename(mask, row, id):
   
   return '{}_spec1D/{}_{}_{}.fits'.format(mask, mask, row, id)
   
def getspectrum(mask, row, id):
   
   filename = getspecfilename(mask, row, id)

   spec = Table.read(filename)
   return spec

def gaussian(x, amp, cen, wid, a, b):
    """1-d gaussian: gaussian(x, amp, cen, wid)"""
    return (amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2)) + a*x + b

# Set up the command line argument parser
parser = argparse.ArgumentParser(description='Correct wavelength and redshifts for air-vacuum issue.')
parser.add_argument('-m', metavar='mask name', type=str, help='name of the mask to be summarized', required=True)
args = parser.parse_args()

directory = './{}_spec1D/'.format(args.m)

# Read in the objects file
objects = Table.read(directory + '{}_objects.fits'.format(args.m))


# Remove the alignment box objects
objects = objects[objects['alignbox'] == False]
print(objects)

# Read in the first object
object = objects[0]
spec = getspectrum(args.m, object['row'], object['id'])

waveArray = spec['wave']
errorArray = np.zeros((len(waveArray), len(objects)))
errorArray[:, :] = np.nan
print(errorArray)
print(errorArray.shape)


for i in range(len(objects)):
   print(i)
   object = objects[i]
   try:
      spec = getspectrum(args.m, object['row'], object['id'])
      errorArray[:, i] = spec['error']
   except:
      print('Spectrum not found.')
   
print(errorArray)
error = np.nanmedian(errorArray, axis=1)
print(error.shape)

print(error)

wavecheck = Table()
wavecheck['wave'] = waveArray
wavecheck['medianerror'] = error
wavecheck['sky'] = error**2
wavecheck = wavecheck[(wavecheck['wave'] > 5578.5 - 25) & (wavecheck['wave'] < 5578.5 + 25)]


gmodel = Model(gaussian)
result = gmodel.fit(wavecheck['sky'], x=wavecheck['wave'], amp=np.max(wavecheck['sky']), cen=5578.5, wid=5, a=0, b=0)
print(result.fit_report())

mu = result.best_values['cen']
FWHM_kms = result.best_values['wid']*2.355/mu*300e3
FWHM_Ang = result.best_values['wid']*2.355

fig, ax = plt.subplots(1, figsize=(10, 7))


ax.plot(wavecheck['wave'], wavecheck['sky'], color='black', drawstyle='steps-mid', label=r'$\rm median\ error$')
ax.plot(wavecheck['wave'], result.best_fit, color='orange', label=r'$\mu = {:0.2f}, \rm FWHM={:0.2f} {}, {:0.0f} {}$'.format(mu, FWHM_Ang, '\ \AA', FWHM_kms, '\ km/s'))
ax.axvline(5578.5, color='blue', label=r'$\rm 5578.5\ (vacuum)$')
ax.axvline(5576.95, color='red', linestyle=':', label=r'$5\rm 5576.95\ (air)$')

ax.legend()
ax.minorticks_on()
ax.set_xlabel(r'$\rm observed\ wavelength\ [\AA]$')
ax.set_ylabel(r'$\rm error\ array\ [arbitrary\ unit]$')
fig.tight_layout()
plt.savefig('{}_wavecheck.pdf'.format(args.m))
plt.close()