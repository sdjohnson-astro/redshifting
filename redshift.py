import os
import numpy as np
#from astropy.table import Table, Column, Row
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline, UnivariateSpline
from lmfit import Model, Parameters
import corner
import time
import warnings
from matplotlib import pyplot as plt

# Set constants
c_kms = 299792.458

# Read in the galaxy eigenspectra
eigen_galaxy = np.load(os.environ['REDSHIFTING'] + '/eigenspectra/eigen_galaxy.npy')

# Interpolate the Eigenspectra for galaxies
eigen_galaxy1_interp = interp1d(eigen_galaxy['wave'], eigen_galaxy['flux1'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_galaxy2_interp = interp1d(eigen_galaxy['wave'], eigen_galaxy['flux2'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_galaxy3_interp = interp1d(eigen_galaxy['wave'], eigen_galaxy['flux3'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_galaxy4_interp = interp1d(eigen_galaxy['wave'], eigen_galaxy['flux4'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)


# Read in the quasar eigenspectra
eigen_qso = np.load(os.environ['REDSHIFTING'] + '/eigenspectra/eigen_qso.npy')
eigen_qso1_interp = interp1d(eigen_qso['wave'], eigen_qso['flux1'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_qso2_interp = interp1d(eigen_qso['wave'], eigen_qso['flux2'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_qso3_interp = interp1d(eigen_qso['wave'], eigen_qso['flux3'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_qso4_interp = interp1d(eigen_qso['wave'], eigen_qso['flux4'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)


# Read in the star eigenspectra
eigen_star = np.load(os.environ['REDSHIFTING'] + '/eigenspectra/eigen_star.npy')
eigen_star1_interp = interp1d(eigen_star['wave'], eigen_star['flux1'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_star2_interp = interp1d(eigen_star['wave'], eigen_star['flux2'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_star3_interp = interp1d(eigen_star['wave'], eigen_star['flux3'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_star4_interp = interp1d(eigen_star['wave'], eigen_star['flux4'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_star5_interp = interp1d(eigen_star['wave'], eigen_star['flux5'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_star6_interp = interp1d(eigen_star['wave'], eigen_star['flux6'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_star7_interp = interp1d(eigen_star['wave'], eigen_star['flux7'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_star8_interp = interp1d(eigen_star['wave'], eigen_star['flux8'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_star9_interp = interp1d(eigen_star['wave'], eigen_star['flux9'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_star10_interp = interp1d(eigen_star['wave'], eigen_star['flux10'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_star11_interp = interp1d(eigen_star['wave'], eigen_star['flux11'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)


# Read in quasar HW templates
quasar_HW = np.load(os.environ['REDSHIFTING'] + '/eigenspectra/quasar_HW.npy')
quasar_HW1_interp = interp1d(quasar_HW['wave'], quasar_HW['flux1'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
quasar_HW2_interp = interp1d(quasar_HW['wave'], quasar_HW['flux2'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)



# Read in LATIS z=2-3 galaxy 
eigen_latis = np.load(os.environ['REDSHIFTING'] + '/eigenspectra/eigen_latis.npy')
eigen_latis1_interp = interp1d(eigen_latis['wave'], eigen_latis['flux1'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_latis2_interp = interp1d(eigen_latis['wave'], eigen_latis['flux2'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_latis3_interp = interp1d(eigen_latis['wave'], eigen_latis['flux3'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_latis4_interp = interp1d(eigen_latis['wave'], eigen_latis['flux4'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)
eigen_latis5_interp = interp1d(eigen_latis['wave'], eigen_latis['flux5'], fill_value=(-999.0, -999.0), assume_sorted=False, bounds_error=False)


# calculat a redshift plus a velocity in km/s
def dz(z, dv):
   znew = dv/c_kms*(1.0 + z) + z

   return znew

# calculate the velocity difference in km/s between two redshfits
def dv(z0, z1):
   
   dv = (z1 - z0)/(1.0 + z0)*c_kms
   
   return dv


# Return a galaxy spectrum model evaluated on wavelength grid wave
# at redshift z with SDSS eigenspectra coefficients eigen1, eigen2, eigen3, eigen4
# Optional low-order polynomials have coffieicnts fluxcal0, fluxcal1, fluxcal2
# where 0, 1, and 2 are for 0th, 1st, and 2nd order polynomials
# There are better ways to implement this last option
def eigensum_galaxy(wave, z, eigen1, eigen2, eigen3, eigen4,
                    fluxcal0, fluxcal1, fluxcal2):
   
   wave_zp1 = wave/(1.0 + z)
   flux = eigen_galaxy1_interp(wave_zp1)*eigen1 \
          + eigen_galaxy2_interp(wave_zp1)*eigen2 \
          + eigen_galaxy3_interp(wave_zp1)*eigen3 \
          + eigen_galaxy4_interp(wave_zp1)*eigen4 \
          + fluxcal0 + fluxcal1*wave + fluxcal2*wave**2
   
   return flux
   

# Same as eigensum_galaxy but with the SDSS quasar eigenspectra
def eigensum_qso(wave, z, eigen1, eigen2, eigen3, eigen4,
                    fluxcal0, fluxcal1, fluxcal2):
   
   wave_zp1 = wave/(1.0 + z)
   flux = eigen_qso1_interp(wave_zp1)*eigen1 \
          + eigen_qso2_interp(wave_zp1)*eigen2 \
          + eigen_qso3_interp(wave_zp1)*eigen3 \
          + eigen_qso4_interp(wave_zp1)*eigen4 \
          + fluxcal0 + fluxcal1*wave + fluxcal2*wave**2
   
   return flux
   
# Same as eigensum_galaxy but with the LATIS galaxy templates
def eigensum_latis(wave, z, eigen1, eigen2, eigen3, eigen4, eigen5, 
                    fluxcal0, fluxcal1, fluxcal2):
   
   wave_zp1 = wave/(1.0 + z)
   flux =   eigen_latis1_interp(wave_zp1)*eigen1 \
          + eigen_latis2_interp(wave_zp1)*eigen2 \
          + eigen_latis3_interp(wave_zp1)*eigen3 \
          + eigen_latis4_interp(wave_zp1)*eigen4 \
          + eigen_latis5_interp(wave_zp1)*eigen5 \
          + fluxcal0 + fluxcal1*wave + fluxcal2*wave**2
   
   return flux

# Same as eigensum_galaxy but with SDSS stellar spectral templates.
def eigensum_star(wave, z, eigen1, eigen2, eigen3, eigen4, eigen5, eigen6, eigen7,
                 eigen8, eigen9, eigen10, eigen11,
                 fluxcal0, fluxcal1, fluxcal2):
   
   wave_zp1 = wave/(1.0 + z)
   flux = eigen_star1_interp(wave_zp1)*eigen1 \
          + eigen_star2_interp(wave_zp1)*eigen2 \
          + eigen_star3_interp(wave_zp1)*eigen3 \
          + eigen_star4_interp(wave_zp1)*eigen4 \
          + eigen_star5_interp(wave_zp1)*eigen5 \
          + eigen_star6_interp(wave_zp1)*eigen6 \
          + eigen_star7_interp(wave_zp1)*eigen7 \
          + eigen_star8_interp(wave_zp1)*eigen8 \
          + eigen_star9_interp(wave_zp1)*eigen9 \
          + eigen_star10_interp(wave_zp1)*eigen10 \
          + eigen_star11_interp(wave_zp1)*eigen11 \
          + fluxcal0 + fluxcal1*wave + fluxcal2*wave**2
   
   return flux

# Get best-fit galaxy model at a redshift. flux calibration polynomials not currently included
def fitatz_galaxy(spec, z):
   
   
   # Evaluate eigenspecra at needed redshifts
   
   constant = np.ones(len(spec))
   linear = np.arange(len(spec))
   square = np.arange(len(spec))**2
   
   wave_zp1 = spec['wave']/(1.0 + z)
   eigen1 = eigen_galaxy1_interp(wave_zp1)
   eigen2 = eigen_galaxy2_interp(wave_zp1)
   eigen3 = eigen_galaxy3_interp(wave_zp1)
   eigen4 = eigen_galaxy4_interp(wave_zp1)
   
   mask = (eigen1 != -999.0).astype(float)*spec['mask']
   
   At = np.matrix([eigen1, eigen2, eigen3, eigen4, constant, linear, square])
   A = At.transpose()
   
   one_over_sigmasquared = 1/spec['error']**2
   one_over_sigmasquared[~np.isfinite(one_over_sigmasquared)] = 0.0
   Ci = np.diag(one_over_sigmasquared)
   A = At.transpose()
   Y = np.matrix(spec['flux']).transpose()
   


   AtCi = np.matrix([eigen1*one_over_sigmasquared*mask,
                     eigen2*one_over_sigmasquared*mask,
                     eigen3*one_over_sigmasquared*mask,
                     eigen4*one_over_sigmasquared*mask,
                     constant*one_over_sigmasquared*mask,
                     linear*one_over_sigmasquared*mask,
                     square*one_over_sigmasquared*mask])#At*Ci
   eigenvalues = np.linalg.inv(AtCi*A)*AtCi*Y
   eigenvalues = eigenvalues.getA1()

   model = eigenvalues[0]*eigen1 + eigenvalues[1]*eigen2 \
           + eigenvalues[2]*eigen3 + eigenvalues[3]*eigen4 \
           + eigenvalues[4]*constant + eigenvalues[5]*linear + eigenvalues[6]*square
        
   chi2 = np.sum(np.square(spec['flux'] - model)*one_over_sigmasquared*mask)
   dof = np.sum(mask) - 7
   chi2_pdf = chi2/dof

   return (eigenvalues, model, chi2_pdf)
   
   
   
def fitatz_latis(spec, z):
   
   
   # Evaluate eigenspecra at needed redshifts
   
   constant = np.ones(len(spec))
   linear = np.arange(len(spec))
   square = np.arange(len(spec))**2
   
   wave_zp1 = spec['wave']/(1.0 + z)
   eigen1 = eigen_latis1_interp(wave_zp1)
   eigen2 = eigen_latis2_interp(wave_zp1)
   eigen3 = eigen_latis3_interp(wave_zp1)
   eigen4 = eigen_latis4_interp(wave_zp1)
   eigen5 = eigen_latis5_interp(wave_zp1)
   
   mask = (eigen1 != -999.0).astype(float)*spec['mask']
   
   At = np.matrix([eigen1, eigen2, eigen3, eigen4, eigen5, constant, linear, square])
   A = At.transpose()
   
   one_over_sigmasquared = 1/spec['error']**2
   one_over_sigmasquared[~np.isfinite(one_over_sigmasquared)] = 0.0
   Ci = np.diag(one_over_sigmasquared)
   A = At.transpose()
   Y = np.matrix(spec['flux']).transpose()
   


   AtCi = np.matrix([eigen1*one_over_sigmasquared*mask,
                     eigen2*one_over_sigmasquared*mask,
                     eigen3*one_over_sigmasquared*mask,
                     eigen4*one_over_sigmasquared*mask,
                     eigen5*one_over_sigmasquared*mask,
                     constant*one_over_sigmasquared*mask,
                     linear*one_over_sigmasquared*mask,
                     square*one_over_sigmasquared*mask])#At*Ci
   eigenvalues = np.linalg.inv(AtCi*A)*AtCi*Y
   eigenvalues = eigenvalues.getA1()

   model = eigenvalues[0]*eigen1 + eigenvalues[1]*eigen2 \
           + eigenvalues[2]*eigen3 + eigenvalues[3]*eigen4 + eigenvalues[4]*eigen5 \
           + eigenvalues[5]*constant + eigenvalues[6]*linear + eigenvalues[7]*square
   model = model*mask     
   chi2 = np.sum(np.square(spec['flux'] - model)*one_over_sigmasquared*mask)
   dof = np.sum(mask) - 8
   chi2_pdf = chi2/dof

   return (eigenvalues, model, chi2_pdf)
   
   
      
 # Same as fitatz_galaxy but for quasars     
def fitatz_qso(spec, z):
   
      constant = np.ones(len(spec))
      linear = np.arange(len(spec))
      square = np.arange(len(spec))**2
   
      # Evaluate eigenspecra at needed redshifts
      wave_zp1 = spec['wave']/(1.0 + z)
      eigen1 = eigen_qso1_interp(wave_zp1)
      eigen2 = eigen_qso2_interp(wave_zp1)
      eigen3 = eigen_qso3_interp(wave_zp1)
      eigen4 = eigen_qso4_interp(wave_zp1)
      
      mask = (eigen1 != -999.0).astype(float)
      
      At = np.matrix([eigen1, eigen2, eigen3, eigen4, constant, linear, square])
      C = np.diag(spec['error']**2)
      one_over_sigmasquared = 1/spec['error']**2
      one_over_sigmasquared[~np.isfinite(one_over_sigmasquared)] = 0.0
      
      index = np.where((spec['mask'] == 0) | (spec['error'] == 0.0))
      one_over_sigmasquared[index] = 0.0
      
      Ci = np.diag(1/spec['error']**2)
      A = At.transpose()
      Y = np.matrix(spec['flux']).transpose()
      
      
      
      AtCi = np.matrix([eigen1*one_over_sigmasquared*mask,       
                        eigen2*one_over_sigmasquared*mask,
                        eigen3*one_over_sigmasquared*mask,
                        eigen4*one_over_sigmasquared*mask,
                        constant*one_over_sigmasquared*mask,
                        linear*one_over_sigmasquared*mask,
                        square*one_over_sigmasquared*mask])#At*Ci
      eigenvalues = np.linalg.inv(AtCi*A)*AtCi*Y
      eigenvalues = eigenvalues.getA1()
      
      model = eigenvalues[0]*eigen1 + eigenvalues[1]*eigen2 \
              + eigenvalues[2]*eigen3 + eigenvalues[3]*eigen4 \
              + eigenvalues[4]*constant + eigenvalues[5]*linear + eigenvalues[6]*square
              
      chi2 = np.sum(np.square(spec['flux'] - model)*one_over_sigmasquared*mask)
      dof = np.sum(mask) - 7.0
      chi2_pdf = chi2/dof
      
      return (eigenvalues, model, chi2_pdf)
      
 # Same as fitatz_galaxy but for stars     
def fitatz_star(spec, z):
   
      constant = np.ones(len(spec))
      linear = np.arange(len(spec))
      square = np.arange(len(spec))**2
   
      # Precompute some matrices
      C = np.diag(spec['error']**2)
      Ci = np.diag(1/spec['error']**2)
      one_over_sigmasquared = 1/spec['error']**2
      one_over_sigmasquared[~np.isfinite(one_over_sigmasquared)] = 0.0
      
      index = np.where((spec['mask'] == 0) | (spec['error'] == 0.0))
      one_over_sigmasquared[index] = 0.0
      
      wave = spec['wave']
      flux = spec['flux']
      error = spec['error']
      Y = np.matrix(flux).transpose()
      
      # Evaluate eigenspecra at needed redshifts
      wave_zp1 = wave/(1.0 + z)
      eigen1 = eigen_star1_interp(wave_zp1)
      eigen2 = eigen_star2_interp(wave_zp1)
      eigen3 = eigen_star3_interp(wave_zp1)
      eigen4 = eigen_star4_interp(wave_zp1)
      eigen5 = eigen_star5_interp(wave_zp1)
      eigen6 = eigen_star6_interp(wave_zp1)
      eigen7 = eigen_star7_interp(wave_zp1)
      eigen8 = eigen_star8_interp(wave_zp1)
      eigen9 = eigen_star9_interp(wave_zp1)
      eigen10 = eigen_star10_interp(wave_zp1)
      eigen11 = eigen_star11_interp(wave_zp1)
      
      mask = (eigen1 != -999.0).astype(float)
      
      
   
      At = np.matrix([eigen1, eigen2, eigen3, eigen4, eigen5, eigen6, eigen7,
                      eigen8, eigen9, eigen10, eigen11, constant, linear, square])
      A = At.transpose()
      
   
   
      AtCi = np.matrix([eigen1*one_over_sigmasquared*mask,
                        eigen2*one_over_sigmasquared*mask,
                        eigen3*one_over_sigmasquared*mask,
                        eigen4*one_over_sigmasquared*mask,
                        eigen5*one_over_sigmasquared*mask,
                        eigen6*one_over_sigmasquared*mask,
                        eigen7*one_over_sigmasquared*mask,
                        eigen8*one_over_sigmasquared*mask,
                        eigen9*one_over_sigmasquared*mask,
                        eigen10*one_over_sigmasquared*mask,
                        eigen11*one_over_sigmasquared*mask,
                        constant*one_over_sigmasquared*mask,
                        linear*one_over_sigmasquared*mask,
                        square*one_over_sigmasquared*mask])#At*Ci
      eigenvalues = np.linalg.inv(AtCi*A)*AtCi*Y
      eigenvalues = eigenvalues.getA1()
   
      model = eigenvalues[0]*eigen1 + eigenvalues[1]*eigen2 \
              + eigenvalues[2]*eigen3 + eigenvalues[3]*eigen4 \
              + eigenvalues[4]*eigen5 + eigenvalues[5]*eigen6 \
              + eigenvalues[6]*eigen7 + eigenvalues[7]*eigen8 \
              + eigenvalues[8]*eigen9 + eigenvalues[9]*eigen10 \
              + eigenvalues[10]*eigen11 \
              + eigenvalues[11]*constant + eigenvalues[12]*linear + eigenvalues[13]*square
              
      model = model*mask     
      chi2 = np.sum(np.square(flux - model)*one_over_sigmasquared*mask)
      
      dof = np.sum(mask) - 13.0
      chi2_pdf = chi2/dof
      
      
      return (eigenvalues, model, chi2_pdf)
      
      
 # Same as fitatz_galaxy but for quasars with Hewitt and Wild cross-correlation template and flux calibration polynomials
def fitatz_qso_hw_poly(spec, z):
   
   
   C = np.diag(spec['error']**2)
   Ci = np.diag(1/spec['error']**2)
   one_over_sigmasquared = 1/spec['error']**2
   
   index = np.where((spec['mask'] == 0) | (spec['error'] == 0.0))
   one_over_sigmasquared[index] = 0.0
   
   wave = spec['wave']
   flux = spec['flux']
   error = spec['error']
   Y = np.matrix(flux).transpose()
   
   
   constant = np.ones(len(spec))
   linear = np.arange(len(spec))
   square = np.arange(len(spec))**2
   
   
   wave_zp1 = wave/(1.0 + z)
   eigen1 = quasar_HW1_interp(wave_zp1)
   eigen2 = quasar_HW2_interp(wave_zp1)
   mask = (eigen1 != -999.0).astype(float)*spec['mask']


   At = np.matrix([eigen1, eigen2, constant, linear, square])
   A = At.transpose()
   


   AtCi = np.matrix([eigen1*one_over_sigmasquared,
                     eigen2*one_over_sigmasquared,
                     constant*one_over_sigmasquared,
                     linear*one_over_sigmasquared,
                     square*one_over_sigmasquared])#At*Ci
   eigenvalues = np.linalg.inv(AtCi*A)*AtCi*Y
   eigenvalues = eigenvalues.getA1()

   model_qso = eigenvalues[0]*eigen1 + eigenvalues[1]*eigen2
   model_poly = eigenvalues[2]*constant + eigenvalues[3]*linear + eigenvalues[4]*square
   model = model_qso + model_poly
                 
   chi2 = np.sum(np.square(flux - model)*one_over_sigmasquared)
   dof = np.sum(mask) - 5.0
   chi2_pdf = chi2/dof
   
   return (eigenvalues, model, chi2_pdf)
   
      

# Return chi2, chi2pdf, and best-fit coefficients for stellar model fit to spectrum over a grid of wavelengths from zim to zmax with steps of dz
# The spectrum *must* my a numpy nparray or equivalent
# with wave, flux, error, and mask columns
# wavelength in angstroms
# flux and error in flambda
# mask = 1 for good data, = 0 for bad data               
def findz_star(spec, zmin=-0.005, zmax=0.005, dz=0.000025):
   
   zs = np.arange(zmin, zmax, dz)
   zeros = np.zeros(len(zs))
   indices = np.arange(0, len(zs), dtype=int)
   redshifts = np.ndarray(len(zs),
                         dtype={'names':('index', 'z', 'chi2', 'chi2_pdf', 'eigen1', 'eigen2',
                                'eigen3', 'eigen4', 'eigen5', 'eigen6', 'eigen7',
                                'eigen8', 'eigen9', 'eigen10', 'eigen11',
                                'fluxcal0', 'fluxcal1', 'fluxcal2'), 
                                'formats':(int, float, float, float, float, float,
                                           float, float, float, float, float,
                                           float, float, float, float,
                                           float, float, float)})
   redshifts['index'] = indices
   redshifts['z'] = zs                                                         
   redshifts['eigen1']   = 0.0   
   redshifts['eigen2']   = 0.0
   redshifts['eigen3']   = 0.0
   redshifts['eigen4']   = 0.0
   redshifts['eigen5']   = 0.0
   redshifts['eigen6']   = 0.0
   redshifts['eigen7']   = 0.0
   redshifts['eigen8']   = 0.0
   redshifts['eigen9']   = 0.0
   redshifts['eigen10']   = 0.0
   redshifts['eigen11']   = 0.0
   redshifts['fluxcal0'] = 0.0
   redshifts['fluxcal1'] = 0.0
   redshifts['fluxcal2'] = 0.0
   
   constant = np.ones(len(spec))
   linear = np.arange(len(spec))
   square = np.arange(len(spec))**2
   
   #redshifts = Table(redshifts)
   
   # Precompute some matrices
   C = np.diag(spec['error']**2)
   Ci = np.diag(1/spec['error']**2)
   one_over_sigmasquared = 1/spec['error']**2
   
   index = np.where((spec['mask'] == 0) | (spec['error'] == 0.0))
   one_over_sigmasquared[index] = 0.0
   
   wave = spec['wave']
   flux = spec['flux']
   error = spec['error']
   Y = np.matrix(flux).transpose()
   
   def fitatz_star_precompute(z):
   
         # Evaluate eigenspecra at needed redshifts
         wave_zp1 = wave/(1.0 + z)
         eigen1 = eigen_star1_interp(wave_zp1)
         eigen2 = eigen_star2_interp(wave_zp1)
         eigen3 = eigen_star3_interp(wave_zp1)
         eigen4 = eigen_star4_interp(wave_zp1)
         eigen5 = eigen_star5_interp(wave_zp1)
         eigen6 = eigen_star6_interp(wave_zp1)
         eigen7 = eigen_star7_interp(wave_zp1)
         eigen8 = eigen_star8_interp(wave_zp1)
         eigen9 = eigen_star9_interp(wave_zp1)
         eigen10 = eigen_star10_interp(wave_zp1)
         eigen11 = eigen_star11_interp(wave_zp1)
         
         mask = (eigen1 != -999.0).astype(float)
      
         At = np.matrix([eigen1, eigen2, eigen3, eigen4, eigen5, eigen6, eigen7,
                         eigen8, eigen9, eigen10, eigen11, constant, linear, square])
         A = At.transpose()
         
      
      
         AtCi = np.matrix([eigen1*one_over_sigmasquared*mask,
                           eigen2*one_over_sigmasquared*mask,
                           eigen3*one_over_sigmasquared*mask,
                           eigen4*one_over_sigmasquared*mask,
                           eigen5*one_over_sigmasquared*mask,
                           eigen6*one_over_sigmasquared*mask,
                           eigen7*one_over_sigmasquared*mask,
                           eigen8*one_over_sigmasquared*mask,
                           eigen9*one_over_sigmasquared*mask,
                           eigen10*one_over_sigmasquared*mask,
                           eigen11*one_over_sigmasquared*mask,
                           constant*one_over_sigmasquared*mask,
                           linear*one_over_sigmasquared*mask,
                           square*one_over_sigmasquared*mask])#At*Ci
         eigenvalues = np.linalg.inv(AtCi*A)*AtCi*Y
         eigenvalues = eigenvalues.getA1()
      
         model = eigenvalues[0]*eigen1 + eigenvalues[1]*eigen2 \
                 + eigenvalues[2]*eigen3 + eigenvalues[3]*eigen4 \
                 + eigenvalues[4]*eigen5 + eigenvalues[5]*eigen6 \
                 + eigenvalues[6]*eigen7 + eigenvalues[7]*eigen8 \
                 + eigenvalues[8]*eigen9 + eigenvalues[9]*eigen10 \
                 + eigenvalues[10]*eigen11 \
                 + + eigenvalues[11]*constant + eigenvalues[12]*linear + eigenvalues[13]*square
         model = model*mask
              
         chi2 = np.sum(np.square(flux - model)*one_over_sigmasquared*mask)
         dof = np.sum(mask)*1.0 - 11.0
         chi2_pdf = chi2/dof
         return (eigenvalues, chi2, chi2_pdf)

   for redshift in redshifts:
      
      
      eigenvalues, chi2, chi2_pdf = fitatz_star_precompute(redshift['z'])
      redshift['chi2'] = chi2
      redshift['chi2_pdf'] = chi2_pdf
      redshift['eigen1'] = eigenvalues[0]
      redshift['eigen2'] = eigenvalues[1]
      redshift['eigen3'] = eigenvalues[2]
      redshift['eigen4'] = eigenvalues[3]
      redshift['eigen5'] = eigenvalues[4]
      redshift['eigen6'] = eigenvalues[5]
      redshift['eigen7'] = eigenvalues[6]
      redshift['eigen8'] = eigenvalues[7]
      redshift['eigen9'] = eigenvalues[8]
      redshift['eigen10'] = eigenvalues[9]
      redshift['eigen11'] = eigenvalues[10]
      redshifts['fluxcal0'] = eigenvalues[11]
      redshifts['fluxcal1'] = eigenvalues[12]
      redshifts['fluxcal2'] = eigenvalues[13]
   
   return redshifts
   
   
# Return chi2, chi2pdf, and best-fit coefficients for quasar model fit to spectrum over a grid of wavelengths from zim to zmax with steps of dz
# The spectrum *must* my a numpy nparray or equivalent
# with wave, flux, error, and mask columns
# wavelength in angstroms
# flux and error in flambda
# mask = 1 for good data, = 0 for bad data   
# Does *not* include IGM attenuation
def findz_qso(spec, zmin=0.0, zmax=3.0, dz=0.0005):
   
   
   
   # Create an array of Lya absorption
   zs = np.arange(zmin, zmax, dz)
   zeros = np.zeros(len(zs))
   indices = np.arange(0, len(zs), dtype=int)
   redshifts = np.ndarray(len(zs),
                         dtype={'names':('index', 'z', 'chi2', 'chi2_pdf', 'eigen1', 'eigen2',
                                'eigen3', 'eigen4', 'fluxcal0', 'fluxcal1', 'fluxcal2'), 
                                'formats':(int, float, float, float, float, float,
                                           float, float, float, float, float)})
   redshifts['index'] = indices
   redshifts['z'] = zs                                                         
   redshifts['eigen1']   = 0.0   
   redshifts['eigen2']   = 0.0
   redshifts['eigen3']   = 0.0
   redshifts['eigen4']   = 0.0
   redshifts['fluxcal0'] = 0.0
   redshifts['fluxcal1'] = 0.0
   redshifts['fluxcal2'] = 0.0
   
   constant = np.ones(len(spec))
   linear = np.arange(len(spec))
   square = np.arange(len(spec))**2
   
   
   
   # Precompute some matrices
   C = np.diag(spec['error']**2)
   Ci = np.diag(1/spec['error']**2)
   one_over_sigmasquared = 1/spec['error']**2
   
   index = np.where((spec['mask'] == 0) | (spec['error'] == 0.0))
   one_over_sigmasquared[index] = 0.0
   
   wave = spec['wave']
   flux = spec['flux']
   error = spec['error']
   Y = np.matrix(flux).transpose()
   
   def fitatz_qso_precompute(z):
   
         # Evaluate eigenspecra at needed redshifts
         wave_zp1 = wave/(1.0 + z)
         eigen1 = eigen_qso1_interp(wave_zp1)
         eigen2 = eigen_qso2_interp(wave_zp1)
         eigen3 = eigen_qso3_interp(wave_zp1)
         eigen4 = eigen_qso4_interp(wave_zp1)
         
         mask = (eigen1 != -999.0).astype(float)*spec['mask']
      
         At = np.matrix([eigen1, eigen2, eigen3, eigen4, constant, linear, square])
         A = At.transpose()
         
      
      
         AtCi = np.matrix([eigen1*one_over_sigmasquared*mask,
                           eigen2*one_over_sigmasquared*mask,
                           eigen3*one_over_sigmasquared*mask,
                           eigen4*one_over_sigmasquared*mask,
                           constant*one_over_sigmasquared*mask,
                           linear*one_over_sigmasquared*mask,
                           square*one_over_sigmasquared*mask])#At*Ci
         eigenvalues = np.linalg.inv(AtCi*A)*AtCi*Y
         eigenvalues = eigenvalues.getA1()
      
         model = eigenvalues[0]*eigen1 + eigenvalues[1]*eigen2 \
                 + eigenvalues[2]*eigen3 + eigenvalues[3]*eigen4 \
                 + eigenvalues[4]*constant + eigenvalues[5]*linear + eigenvalues[6]*square
         model = model*mask     
         chi2 = np.sum(np.square(flux - model)*one_over_sigmasquared*mask)
         dof = np.sum(mask) - 4.0
         chi2_pdf = chi2/dof
         return (eigenvalues, chi2, chi2_pdf)

   for redshift in redshifts:
      
      
      eigenvalues, chi2, chi2_pdf = fitatz_qso_precompute(redshift['z'])
      redshift['chi2'] = chi2
      redshift['chi2_pdf'] = chi2_pdf
      redshift['eigen1'] = eigenvalues[0]
      redshift['eigen2'] = eigenvalues[1]
      redshift['eigen3'] = eigenvalues[2]
      redshift['eigen4'] = eigenvalues[3]
      redshifts['fluxcal0'] = eigenvalues[4]
      redshifts['fluxcal1'] = eigenvalues[5]
      redshifts['fluxcal2'] = eigenvalues[6]

   
   return redshifts
   
   
   
# Return chi2, chi2pdf, and best-fit coefficients for quasar model fit to spectrum over a grid of wavelengths from zim to zmax with steps of dz
# The spectrum *must* my a numpy nparray or equivalent
# with wave, flux, error, and mask columns
# wavelength in angstroms
# flux and error in flambda
# mask = 1 for good data, = 0 for bad data           
# Does *not* include IGM attenuation
def findz_qso_hw(spec, zmin=0.0, zmax=3.0, dz=0.0005):
   
   
   
   # Create an array of Lya absorption
   zs = np.arange(zmin, zmax, dz)
   zeros = np.zeros(len(zs))
   indices = np.arange(0, len(zs), dtype=int)
   redshifts = np.ndarray(len(zs),
                         dtype={'names':('index', 'z', 'chi2', 'chi2_pdf', 'eigen1', 'eigen2',
                                'fluxcal0', 'fluxcal1', 'fluxcal2'), 
                                'formats':(int, float, float, float, float, float,
                                           float, float, float, float, float)})
   redshifts['index'] = indices
   redshifts['z'] = zs                                                         
   redshifts['eigen1']   = 0.0   
   redshifts['eigen2']   = 0.0
   redshifts['fluxcal0'] = 0.0
   redshifts['fluxcal1'] = 0.0
   redshifts['fluxcal2'] = 0.0
   
   #redshifts = Table(redshifts)
   
   # Precompute some matrices
   C = np.diag(spec['error']**2)
   Ci = np.diag(1/spec['error']**2)
   one_over_sigmasquared = 1/spec['error']**2
   
   index = np.where((spec['mask'] == 0) | (spec['error'] == 0.0))
   one_over_sigmasquared[index] = 0.0
   
   wave = spec['wave']
   flux = spec['flux']
   error = spec['error']
   Y = np.matrix(flux).transpose()
   
   
   constant = np.ones(len(spec))
   linear = np.arange(len(spec))
   square = np.arange(len(spec))**2
   
   def fitatz_qso_hw_precompute(z):
   
         # Evaluate eigenspecra at needed redshifts
         wave_zp1 = wave/(1.0 + z)
         eigen1 = quasar_HW1_interp(wave_zp1)
         eigen2 = quasar_HW2_interp(wave_zp1)
         mask = (eigen1 != -999.0).astype(float)*spec['mask']

      
         At = np.matrix([eigen1, eigen2, constant, linear, square])
         A = At.transpose()
         
      
      
         AtCi = np.matrix([eigen1*one_over_sigmasquared,
                           eigen2*one_over_sigmasquared,
                           constant*one_over_sigmasquared,
                           linear*one_over_sigmasquared,
                           square*one_over_sigmasquared])#At*Ci
         eigenvalues = np.linalg.inv(AtCi*A)*AtCi*Y
         eigenvalues = eigenvalues.getA1()
      
         model_qso = eigenvalues[0]*eigen1 + eigenvalues[1]*eigen2
         model_poly = eigenvalues[2]*constant + eigenvalues[3]*linear + eigenvalues[4]*square
         model = model_qso + model_poly
                       
         chi2 = np.sum(np.square(flux - model)*one_over_sigmasquared)
         dof = np.sum(mask) - 5.0
         chi2_pdf = chi2/dof
         return (eigenvalues, chi2, chi2_pdf)

   for redshift in redshifts:
      
      
      eigenvalues, chi2, chi2_pdf = fitatz_qso_hw_precompute(redshift['z'])
      redshift['chi2'] = chi2
      
      redshift['chi2_pdf'] = chi2_pdf
      
      redshift['eigen1'] = eigenvalues[0]
      redshift['eigen2'] = eigenvalues[1]
      redshift['fluxcal0'] = eigenvalues[2]
      redshift['fluxcal1'] = eigenvalues[3]
      redshift['fluxcal2'] = eigenvalues[4]
   
   return redshifts
   


   
   
   
   
#      # Evaluate eigenspecra at needed redshifts
#      wave_zp1 = spec['wave']/(1.0 + z)
#      hw1_qso = quasar_HW1_interp(wave_zp1)
#      hw2_qso = quasar_HW2_interp(wave_zp1)
#      
#      constant = np.ones(len(spec))
#      linear = np.arange(len(spec))
#      square = np.arange(len(spec))**2
#      
#      At = np.matrix([hw1_qso, hw2_qso,
#                      constant, linear, square])
#      C = np.diag(spec['error']**2)
#      one_over_sigmasquared = 1/spec['error']**2
#      Ci = np.diag(1/spec['error']**2)
#      A = At.transpose()
#      Y = np.matrix(spec['flux']).transpose()
#      
#      
#      
#      AtCi = np.matrix([hw1_qso*one_over_sigmasquared,
#                        hw2_qso*one_over_sigmasquared,
#                        constant*one_over_sigmasquared,
#                        linear*one_over_sigmasquared,
#                        square*one_over_sigmasquared])#At*Ci
#      eigenvalues = np.linalg.inv(AtCi*A)*AtCi*Y
#      eigenvalues = eigenvalues.getA1()
#      
#      model_qso = eigenvalues[0]*hw1_qso + eigenvalues[1]*hw2_qso
#      model_poly = eigenvalues[2]*constant + eigenvalues[3]*linear + eigenvalues[4]*square
#      
#      model = model_qso + model_poly
#              
#      chi2 = np.sum(np.square(spec['flux'] - model)*one_over_sigmasquared)
#      
#      
#      mask = (hw1_qso != -999.0).astype(float)*spec['mask']      
#      
#      dof = np.sum(mask) - 13.0
#      chi2_pdf = chi2/dof
#      
#      return (eigenvalues, model, chi2_pdf)
  
   
# Return chi2, chi2pdf, and best-fit coefficients for stellar model fit to spectrum over a grid of wavelengths from zim to zmax with steps of dz
# The spectrum *must* my a numpy nparray or equivalent
# with wave, flux, error, and mask columns
# wavelength in angstroms
# flux and error in flambda
# mask = 1 for good data, = 0 for bad data           
# Appropriate for low-z galaxies observed in rest-frame optical
# Does not include IGM attenuation
def findz_galaxy(spec, zmin=-0.1, zmax=1.5, dz=0.0001):
   
     
   # Create an array of Lya absorption
   zs = np.arange(zmin, zmax, dz)
   zeros = np.zeros(len(zs))
   indices = np.arange(0, len(zs), dtype=int)
   redshifts = np.ndarray(len(zs),
                         dtype={'names':('index', 'z', 'chi2', 'chi2_pdf', 'eigen1', 'eigen2',
                                'eigen3', 'eigen4', 'fluxcal0', 'fluxcal1', 'fluxcal2'), 
                                'formats':(int, float, float, float, float, float,
                                           float, float, float, float, float)})
   redshifts['index'] = indices
   redshifts['z'] = zs                                                         
   redshifts['eigen1']   = 0.0   
   redshifts['eigen2']   = 0.0
   redshifts['eigen3']   = 0.0
   redshifts['eigen4']   = 0.0
   redshifts['fluxcal0'] = 0.0
   redshifts['fluxcal1'] = 0.0
   redshifts['fluxcal2'] = 0.0
   
   
   constant = np.ones(len(spec))
   linear = np.arange(len(spec))
   square = np.arange(len(spec))**2
   
   #redshifts = Table(redshifts)
   
   # Precompute some matrices
   C = np.diag(spec['error']**2)
   Ci = np.diag(1/spec['error']**2)
   one_over_sigmasquared = 1/spec['error']**2
   
   index = np.where((np.isnan(spec['flux'])) | (np.isnan(spec['error'])))
   spec[index]['flux'] = 0.0
   spec[index]['error'] = 0.0
   spec[index]['mask'] = 0
   
   index = np.where((spec['mask'] == 0) | (spec['error'] == 0.0))
   spec[index]['mask'] = 0
   one_over_sigmasquared[index] = 0.0
   
   index = np.where(spec['mask'] == 0)
   one_over_sigmasquared[index] = 0.0
   
   wave = spec['wave']
   flux = spec['flux']
   error = spec['error']
   Y = np.matrix(flux).transpose()
   
   def fitatz_galaxy_precompute(z):
   
         # Evaluate eigenspecra at needed redshifts
         wave_zp1 = wave/(1.0 + z)
         eigen1 = eigen_galaxy1_interp(wave_zp1)
         eigen2 = eigen_galaxy2_interp(wave_zp1)
         eigen3 = eigen_galaxy3_interp(wave_zp1)
         eigen4 = eigen_galaxy4_interp(wave_zp1)
         
         mask = (eigen1 != -999.0).astype(float)*spec['mask']
         
         At = np.matrix([eigen1, eigen2, eigen3, eigen4, constant, linear, square])
         A = At.transpose()
         
      
      
         AtCi = np.matrix([eigen1*one_over_sigmasquared*mask,
                           eigen2*one_over_sigmasquared*mask,
                           eigen3*one_over_sigmasquared*mask,
                           eigen4*one_over_sigmasquared*mask,
                           constant*one_over_sigmasquared*mask,
                           linear*one_over_sigmasquared*mask,
                           square*one_over_sigmasquared*mask])#At*Ci
         eigenvalues = np.linalg.inv(AtCi*A)*AtCi*Y
         eigenvalues = eigenvalues.getA1()
      
         model = eigenvalues[0]*eigen1 + eigenvalues[1]*eigen2 \
                 + eigenvalues[2]*eigen3 + eigenvalues[3]*eigen4 \
                 + eigenvalues[4]*constant + eigenvalues[5]*linear + eigenvalues[6]*square
              
         chi2 = np.sum(np.square(flux - model)*one_over_sigmasquared*mask)
         dof = np.sum(mask) - 7
         chi2_pdf = chi2/dof
      
         return (eigenvalues, chi2, chi2_pdf)

   for redshift in redshifts:
      
      
      eigenvalues, chi2, chi2_pdf = fitatz_galaxy_precompute(redshift['z'])
      redshift['chi2'] = chi2
      redshift['chi2_pdf'] = chi2_pdf
      redshift['eigen1'] = eigenvalues[0]
      redshift['eigen2'] = eigenvalues[1]
      redshift['eigen3'] = eigenvalues[2]
      redshift['eigen4'] = eigenvalues[3]
      redshifts['fluxcal0'] = eigenvalues[4]
      redshifts['fluxcal1'] = eigenvalues[5]
      redshifts['fluxcal2'] = eigenvalues[6]

   
   return redshifts
   
   

# Return chi2, chi2pdf, and best-fit coefficients for z=2 galaxy model fit to spectrum over a grid of wavelengths from zim to zmax with steps of dz
# The spectrum *must* my a numpy nparray or equivalent
# with wave, flux, error, and mask columns
# wavelength in angstroms
# flux and error in flambda
# mask = 1 for good data, = 0 for bad data           
# Appropriate for low-z galaxies observed in rest-frame optical
# Does not include IGM attenuation?
def findz_latis(spec, zmin=2.5, zmax=5.0, dz=0.0005):
   
   
   
   # Create an array of Lya absorption
   zs = np.arange(zmin, zmax, dz)
   zeros = np.zeros(len(zs))
   indices = np.arange(0, len(zs), dtype=int)
   redshifts = np.ndarray(len(zs),
                         dtype={'names':('index', 'z', 'chi2', 'chi2_pdf', 'eigen1', 'eigen2',
                                'eigen3', 'eigen4', 'eigen5', 'fluxcal0', 'fluxcal1', 'fluxcal2'), 
                                'formats':(int, float, float, float, float, float, float,
                                           float, float, float, float, float)})
   redshifts['index'] = indices
   redshifts['z'] = zs                                                         
   redshifts['eigen1']   = 0.0   
   redshifts['eigen2']   = 0.0
   redshifts['eigen3']   = 0.0
   redshifts['eigen4']   = 0.0
   redshifts['eigen5']   = 0.0
   redshifts['fluxcal0'] = 0.0
   redshifts['fluxcal1'] = 0.0
   redshifts['fluxcal2'] = 0.0
   
   
   constant = np.ones(len(spec))
   linear = np.arange(len(spec))
   square = np.arange(len(spec))**2
   
   #redshifts = Table(redshifts)
   
   # Precompute some matrices
   C = np.diag(spec['error']**2)
   Ci = np.diag(1/spec['error']**2)
   one_over_sigmasquared = 1/spec['error']**2
   
   index = np.where((np.isnan(spec['flux'])) | (np.isnan(spec['error'])))
   spec[index]['flux'] = 0.0
   spec[index]['error'] = 0.0
   spec[index]['mask'] = 0
   
   index = np.where((spec['mask'] == 0) | (spec['error'] == 0.0))
   one_over_sigmasquared[index] = 0.0
   
   index = np.where(spec['mask'] == 0)
   one_over_sigmasquared[index] = 0.0
   
   wave = spec['wave']
   flux = spec['flux']
   error = spec['error']
   Y = np.matrix(flux).transpose()
   
   def fitatz_latis_precompute(z):
   
         # Evaluate eigenspecra at needed redshifts
         wave_zp1 = wave/(1.0 + z)
         eigen1 = eigen_latis1_interp(wave_zp1)
         eigen2 = eigen_latis2_interp(wave_zp1)
         eigen3 = eigen_latis3_interp(wave_zp1)
         eigen4 = eigen_latis4_interp(wave_zp1)
         eigen5 = eigen_latis5_interp(wave_zp1)
         
         mask = (eigen1 != -999.0).astype(float)*spec['mask']
         
         At = np.matrix([eigen1, eigen2, eigen3, eigen4, eigen5, constant, linear, square])
         A = At.transpose()
         
      
      
         AtCi = np.matrix([eigen1*one_over_sigmasquared*mask,
                           eigen2*one_over_sigmasquared*mask,
                           eigen3*one_over_sigmasquared*mask,
                           eigen4*one_over_sigmasquared*mask,
                           eigen5*one_over_sigmasquared*mask,
                           constant*one_over_sigmasquared*mask,
                           linear*one_over_sigmasquared*mask,
                           square*one_over_sigmasquared*mask])#At*Ci
         eigenvalues = np.linalg.inv(AtCi*A)*AtCi*Y
         eigenvalues = eigenvalues.getA1()
      
         model = eigenvalues[0]*eigen1 + eigenvalues[1]*eigen2 \
                 + eigenvalues[2]*eigen3 + eigenvalues[3]*eigen4 + eigenvalues[4]*eigen5 \
                 + eigenvalues[5]*constant + eigenvalues[6]*linear + eigenvalues[7]*square
              
         chi2 = np.sum(np.square(flux - model)*one_over_sigmasquared*mask)
         dof = np.sum(mask) - 8
         chi2_pdf = chi2/dof
      
         return (eigenvalues, chi2, chi2_pdf)

   for redshift in redshifts:
      
      
      eigenvalues, chi2, chi2_pdf = fitatz_latis_precompute(redshift['z'])
      redshift['chi2'] = chi2
      redshift['chi2_pdf'] = chi2_pdf
      redshift['eigen1'] = eigenvalues[0]
      redshift['eigen2'] = eigenvalues[1]
      redshift['eigen3'] = eigenvalues[2]
      redshift['eigen4'] = eigenvalues[3]
      redshift['eigen5'] = eigenvalues[4]
      redshifts['fluxcal0'] = eigenvalues[5]
      redshifts['fluxcal1'] = eigenvalues[6]
      redshifts['fluxcal2'] = eigenvalues[7]

   
   return redshifts



# Fit a galaxy model with a guess redshift and optionally low-order
# polynomials to account for an flux calibration errors
def fitz_galaxy_emcee(spec, zguess, fluxpoly=True, steps=5000, burn=2000, progress=True, printReport=True, saveCorner='', zMin=None, zMax=None):
   
   flux_median = np.median(spec['flux'])
   parameters = Parameters()
   #                    (Name,      Value,            Vary,   Min,  Max,   Expr)
   parameters.add_many(('z',        zguess,           True,  zMin,  zMax,  None),
                       ('eigen1',   flux_median*0.4,  True,  None,  None,  None),
                       ('eigen2',   flux_median*0.4,  True,  None,  None,  None),
                       ('eigen3',   flux_median*0.1,  True,  None,  None,  None),
                       ('eigen4',   flux_median*0.1,  True,  None,  None,  None),
                       ('fluxcal0', 0.0,              fluxpoly,  None,  None,  None),
                       ('fluxcal1', 0.0,              fluxpoly,  None,  None,  None),
                       ('fluxcal2', 0.0,              fluxpoly,  None,  None,  None))
   
   galaxy_model = Model(eigensum_galaxy, missing='drop')
   result = galaxy_model.fit(spec['flux'], wave=spec['wave'], weights=1/spec['error'],
                             params=parameters, missing='drop')
                             
                             
   emcee_kws = dict(steps=steps, burn=burn, is_weighted=True,
                    progress=progress)
   emcee_params = result.params.copy()
   
   result_emcee = galaxy_model.fit(spec['flux'], wave=spec['wave'], weights=1/spec['error'],
                                   params=emcee_params, method='emcee', nan_policy='omit',
                                   missing='drop', fit_kws=emcee_kws, show_titles=True)
   result_emcee.conf_interval   
                                
   # find the maximum likelihood solution
   highest_prob = np.argmax(result_emcee.lnprob)
   hp_loc = np.unravel_index(highest_prob, result_emcee.lnprob.shape)
   mle_soln = result_emcee.chain[hp_loc]
   
   #result_emcee.conf_interval()
   
   if printReport:
      print(result_emcee.fit_report())
      #print(result_emcee.ci_report())
      
      z_marginalized = np.percentile(result_emcee.flatchain['z'], [50])[0]
      zErrUp = np.percentile(result_emcee.flatchain['z'], [84.1])[0] - np.percentile(result_emcee.flatchain['z'], [50])[0]
      zErrDown = np.percentile(result_emcee.flatchain['z'], [50])[0] - np.percentile(result_emcee.flatchain['z'], [15.9])[0]
      print('z = {:0.5f} +{:0.5f} -{:0.5f}'.format(z_marginalized, zErrUp, zErrDown))
      
      interval68 = np.percentile(result_emcee.flatchain['z'], [15.9, 84.1])
      interval95 = np.percentile(result_emcee.flatchain['z'], [2.28, 97.7])
      print('68% C.I:')
      print(interval68)
      
      print('95% C.I:')
      print(interval95)
      
      
      
   
   if saveCorner != '':
      emcee_corner = corner.corner(result_emcee.flatchain, labels=['z', 'eigen1', 'eigen2', 'eigen3', 'eigen4', 'fluxcal0', 'fluxcal1', 'fluxcal2'],
                                   truths=mle_soln)
      emcee_corner.savefig(saveCorner)
      plt.close()
      
                             
   return result_emcee



   
# Fit a galaxy model with a guess redshift and optionally low-order
# polynomials to account for an flux calibration errors
def fitz_galaxy(spec, zguess, fluxpoly=True):
   
   flux_median = np.median(spec['flux'])
   parameters = Parameters()
   #                    (Name,      Value,            Vary,   Min,  Max,   Expr)
   parameters.add_many(('z',        zguess,           True,  None,  None,  None),
                       ('eigen1',   flux_median*0.4,  True,  None,  None,  None),
                       ('eigen2',   flux_median*0.4,  True,  None,  None,  None),
                       ('eigen3',   flux_median*0.1,  True,  None,  None,  None),
                       ('eigen4',   flux_median*0.1,  True,  None,  None,  None),
                       ('fluxcal0', 0.0,              fluxpoly,  None,  None,  None),
                       ('fluxcal1', 0.0,              fluxpoly,  None,  None,  None),
                       ('fluxcal2', 0.0,              fluxpoly,  None,  None,  None))
   
   galaxy_model = Model(eigensum_galaxy, missing='drop')
   result = galaxy_model.fit(spec['flux'], wave=spec['wave'], weights=1/spec['error'],
                             params=parameters, missing='drop')
                             
                             
   emcee_kws = dict(steps=500, burn=200, is_weighted=True,
                    progress=True)
   emcee_params = result.params.copy()
                             
   return result
   


# Fit a quasar model with a guess redshift and optionally low-order
# polynomials to account for an flux calibration errors
def fitz_qso(spec, zguess, fluxpoly=True):
   
   flux_median = np.median(spec['flux'])
   parameters = Parameters()
   #                    (Name,      Value,            Vary,   Min,  Max,   Expr)
   parameters.add_many(('z',        zguess,           True,  None,  None,  None),
                       ('eigen1',   flux_median*0.4,  True,  None,  None,  None),
                       ('eigen2',   flux_median*0.4,  True,  None,  None,  None),
                       ('eigen3',   flux_median*0.1,  True,  None,  None,  None),
                       ('eigen4',   flux_median*0.1,  True,  None,  None,  None),
                       ('fluxcal0', 0.0,              fluxpoly,  None,  None,  None),
                       ('fluxcal1', 0.0,              fluxpoly,  None,  None,  None),
                       ('fluxcal2', 0.0,              fluxpoly,  None,  None,  None))
   
   galaxy_model = Model(eigensum_qso, missing='drop')
   result = galaxy_model.fit(spec['flux'], wave=spec['wave'], weights=1/spec['error'],
                             params=parameters, missing='drop')
   return result
   

# Fit a star model with a guess redshift and optionally low-order
# polynomials to account for an flux calibration errors 
def fitz_star(spec, zguess, fluxpoly=True):
   
   flux_median = np.median(spec['flux'])
   parameters = Parameters()
   #                    (Name,      Value,            Vary,   Min,  Max,   Expr)
   parameters.add_many(('z',        zguess,           True,  None,  None,  None),
                       ('eigen1',   flux_median*0.1,  True,  None,  None,  None),
                       ('eigen2',   flux_median*0.1,  True,  None,  None,  None),
                       ('eigen3',   flux_median*0.1,  True,  None,  None,  None),
                       ('eigen4',   flux_median*0.1,  True,  None,  None,  None),
                       ('eigen5',   flux_median*0.1,  True,  None,  None,  None),
                       ('eigen6',   flux_median*0.1,  True,  None,  None,  None),
                       ('eigen7',   flux_median*0.1,  True,  None,  None,  None),
                       ('eigen8',   flux_median*0.1,  True,  None,  None,  None),
                       ('eigen9',   flux_median*0.1,  True,  None,  None,  None),
                       ('eigen10',   flux_median*0.1,  True,  None,  None,  None),
                       ('eigen11',   flux_median*0.1,  True,  None,  None,  None),
                       ('fluxcal0', 0.0,              fluxpoly,  None,  None,  None),
                       ('fluxcal1', 0.0,              fluxpoly,  None,  None,  None),
                       ('fluxcal2', 0.0,              fluxpoly,  None,  None,  None))
   
   star_model = Model(eigensum_star, missing='drop')
   result = star_model.fit(spec['flux'], wave=spec['wave'], weights=1/spec['error'],
                             params=parameters, missing='drop')
   return result
   
   


   
