import numpy as np
from astropy.io import fits
from astropy.table import Table

def readeigen(filename):
   
   hdulist = fits.open(filename)
   hdu = hdulist[0]
   
   eigen = np.ndarray(hdu.header['NAXIS1'],
                         dtype={'names':('wave', 'loglam', 'chi2',
                                'flux1', 'flux2','flux3', 'flux4'), 
                                'formats':(float, float, float,
                                           float, float, float, float)}) 
   eigen['loglam'] = hdu.header['COEFF0'] + np.arange(0.0, hdu.header['NAXIS1'], 1.0)*hdu.header['COEFF1']
   eigen['wave'] = 10.0**eigen['loglam']
   eigen['flux1'] = hdu.data[0, :]
   eigen['flux2'] = hdu.data[1, :]
   eigen['flux3'] = hdu.data[2, :]
   eigen['flux4'] = hdu.data[3, :]
   
   #eigen = eigen[(eigen['wave'] > 1800) & (eigen['wave'] < 9800)]
   
   hdulist.close()
   
   return eigen
   
   
def readeigenstar(filename):
   
   hdulist = fits.open(filename)
   hdu = hdulist[0]
   
   eigen = np.ndarray(hdu.header['NAXIS1'],
                         dtype={'names':('wave', 'loglam', 'chi2',
                                'flux1', 'flux2','flux3', 'flux4', 'flux5',
                                'flux6', 'flux7', 'flux8', 'flux9', 'flux10', 'flux11'), 
                                'formats':(float, float, float,
                                           float, float, float, float, float,
                                           float, float, float, float, float, float)}) 
   eigen['loglam'] = hdu.header['COEFF0'] + np.arange(0.0, hdu.header['NAXIS1'], 1.0)*hdu.header['COEFF1']
   eigen['wave'] = 10.0**eigen['loglam']
   eigen['flux1'] = hdu.data[0, :]
   eigen['flux2'] = hdu.data[1, :]
   eigen['flux3'] = hdu.data[2, :]
   eigen['flux4'] = hdu.data[3, :]
   eigen['flux5'] = hdu.data[4, :]
   eigen['flux6'] = hdu.data[5, :]
   eigen['flux7'] = hdu.data[6, :]
   eigen['flux8'] = hdu.data[7, :]
   eigen['flux9'] = hdu.data[8, :]
   eigen['flux10'] = hdu.data[9, :]
   eigen['flux11'] = hdu.data[10, :]
   
   hdulist.close()
   
   return eigen
   
   
def readHW(filename):
   
   spec_HW = Table.read(filename)
   eigen = np.ndarray(len(spec_HW),
                         dtype={'names':('wave', 'flux1', 'flux2'), 
                                'formats':(float, float, float)})
   eigen['wave'] = spec_HW['LAMBDA']
   eigen['flux1'] = spec_HW['RFLUX1']
   eigen['flux2'] = spec_HW['RFLUX2']
   
   return eigen
   

def readLATIS(filename_wave, filename_templates):
   
   wave = fits.getdata(filename_wave)
   templates = fits.getdata(filename_templates)
   
   template1 = templates[:, 0]
   template2 = templates[:, 1]
   template3 = templates[:, 2]
   template4 = templates[:, 3]
   template5 = templates[:, 4]

   
   eigen = np.ndarray(len(wave),
                         dtype={'names':('wave', 'flux1', 'flux2', 'flux3', 'flux4', 'flux5'), 
                                'formats':(float, float, float, float, float, float)})
   eigen['wave'] = wave
   eigen['flux1'] = template1
   eigen['flux2'] = template2
   eigen['flux3'] = template3
   eigen['flux4'] = template4
   eigen['flux5'] = template5
   
   return eigen
   

# Read in the galaxy eigenspectrum
eigen_galaxy = readeigen('spEigenGal-56436.fits')
eigen_galaxy['flux2'] = eigen_galaxy['flux2']*-1

# Write as fits and numpy
Table(eigen_galaxy).write('eigen_galaxy.fits', overwrite=True)
np.save('eigen_galaxy.npy', eigen_galaxy)


# Read in the quasar eigenspectrum
eigen_qso = readeigen('spEigenQSO-55732.fits')

# Write as fits and numpy
Table(eigen_qso).write('eigen_qso.fits', overwrite=True)
np.save('eigen_qso.npy', eigen_qso)


# Read in the quasar eigenspectrum
eigen_star = readeigenstar('spEigenStar-51845.fits')

# Write as fits and numpy
Table(eigen_star).write('eigen_star.fits', overwrite=True)
np.save('eigen_star.npy', eigen_star)


# Read in the HW quasar templates
quasar_eigen = readHW('Hewett2010_template.fits')
np.save('quasar_HW.npy', quasar_eigen)

# Read in the LATIS templates for z~2-3 galaxies
latis_eigen = readLATIS('latis_wave.fits', 'latis_templates.fits')
Table(latis_eigen).write('eigen_latis.fits', overwrite=True)
np.save('eigen_latis.npy', latis_eigen)
