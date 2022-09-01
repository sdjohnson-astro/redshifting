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

dv_vacair = 82.68 # km/s at 7000 Ang. At 4000 and 10000 it is 84.76 and 82.46 km/s respectively so constant dv is a very good approximation.
c_kms = 299792.458

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
parser.add_argument('-convert', metavar='convert air to vacuum?', type=bool, help='onvert air to vacuum?', required=False, default=False)
parser.add_argument('-force', metavar='force convert air to vacuum regardless of fit results?', type=bool, help='force convert air to vacuum?', required=False, default=False)

args = parser.parse_args()

directory = './{}_spec1D/'.format(args.m)

# Read in the objects file
objects = Table.read(directory + '{}_objects.fits'.format(args.m))


# Remove the alignment box objects
objects = objects[objects['alignbox'] == False]
print(objects)

# Read in the first object
object = objects[10]
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
wavecheck['medianerror'] = error*1e17
wavecheck['sky'] = (error*1e17)**2
wavecheck = wavecheck[(wavecheck['wave'] > 5578.5 - 25) & (wavecheck['wave'] < 5578.5 + 25)]


gmodel = Model(gaussian)

parameters = Parameters()
parameters.add_many(('amp',     np.max(wavecheck['sky']) - np.min(wavecheck['sky']), True, None,    None,    None),
                    ('cen',     5578.5,                   True, None,  None,  None),
                    ('wid',     3.0,                      True, 0.0,     None,    None),
                    ('a',       0.0,                      True, None,    None,    None),
                    ('b',       np.min(wavecheck['sky']), True, None,    None,    None))

result = gmodel.fit(wavecheck['sky'], params=parameters, x=wavecheck['wave'])
print(result.fit_report())

wave_vac = 5578.5
wave_air = 5576.95

mu = result.best_values['cen']
FWHM_kms = result.best_values['wid']*2.355/mu*300e3
FWHM_Ang = result.best_values['wid']*2.355

fig, ax = plt.subplots(1, figsize=(10, 7))


ax.plot(wavecheck['wave'], wavecheck['sky'], color='black', drawstyle='steps-mid', label=r'$\rm median\ error$')
ax.plot(wavecheck['wave'], result.best_fit, color='orange', label=r'$\mu = {:0.2f}, \rm FWHM={:0.2f} {}, {:0.0f} {}$'.format(mu, FWHM_Ang, '\ \AA', FWHM_kms, '\ km/s'))
ax.axvline(wave_vac, color='blue', label=r'$\rm 5578.5\ (vacuum)$')
ax.axvline(wave_air, color='red', linestyle=':', label=r'$\rm 5576.95\ (air)$')

ax.legend()
ax.minorticks_on()
ax.set_xlabel(r'$\rm observed\ wavelength\ [\AA]$')
ax.set_ylabel(r'$\rm error\ array\ [arbitrary\ unit]$')
fig.tight_layout()
plt.savefig('{}_wavecheck.pdf'.format(args.m))
plt.close()


dW_vacuum = mu - wave_vac
dW_air = mu - wave_air

print('')
print('')

print('{}:'.format(args.m))

inVacuum = np.abs(dW_vacuum) < np.abs(dW_air)

if inVacuum:
   
   print('spec1D wavelengths are in vacuum. dW={:0.3f}'.format(dW_vacuum))

else:
   print('spec1D wavelengths are in air!!!!!! dW={:0.3f}'.format(dW_vacuum))

os.system('open {}_wavecheck.pdf'.format(args.m))


print('FWHM = {:0.2f} Ang or ~{:0.0f} km/s'.format(FWHM_Ang, FWHM_kms))

if ((inVacuum == False) & args.convert) | (args.convert & args.force):
   
   print('')
   print('Converting air wavelength to vacuum and update redshifts')
   
   dirname_spec1D = '{}_spec1D'.format(args.m)
   dirname_spec1D_vac = '{}_spec1D_vac'.format(args.m)
   dirname_spec1D_air = '{}_spec1D_air'.format(args.m)
   
   if (not os.path.exists(dirname_spec1D_vac) & (not os.path.exists(dirname_spec1D_air))):
      
      os.makedirs(dirname_spec1D_vac)
      
      # Read in the objects file
      objects = Table.read('{}_spec1D/{}_objects.fits'.format(args.m, args.m))
      print(objects)
      objects['redshift'] = objects['redshift'] + dv_vacair/c_kms*(1 + objects['redshift'])
      
      objects.write('{}_spec1D_vac/{}_objects.fits'.format(args.m, args.m))
      
      # Loop through all the redshifts files and update.
      filelist_redshifts = glob.glob('{}_spec1D/{}*_redshift.fits'.format(args.m, args.m))
      print(filelist_redshifts)
      for filename_redshifts in filelist_redshifts:
         print(filename_redshifts)
         filename_spec1D = filename_redshifts.replace('_redshift.fits', '.fits')

         print(filename_spec1D)
         print(filename_redshifts)
         
         redshifts = Table.read(filename_redshifts)
         spec = Table.read(filename_spec1D)
         print(redshifts)
         print(spec)
         
         redshifts['z'] = redshifts['z'] + dv_vacair/c_kms*(1 + redshifts['z'])
         redshifts.write(filename_redshifts.replace('_spec1D', '_spec1D_vac'), overwrite=False)
         
         spec['wave'] = airtovac(spec['wave']*u.Angstrom)
         spec.write(filename_spec1D.replace('_spec1D', '_spec1D_vac'), overwrite=False)
         
         
         print('')
      
      os.rename('{}_spec1D'.format(args.m), '{}_spec1D_air'.format(args.m))
      os.rename('{}_spec1D_vac'.format(args.m), '{}_spec1D'.format(args.m))
      
      
   else:
      print('Vacuum directory already exists, will not overwrite')
      

