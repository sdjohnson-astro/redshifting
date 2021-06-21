#This code revises redshifts from spectra that have been remade with better wavecal.
#For Q=2 objects, It applies the masking from the old version and reruns redshift measurement.
# If different is greater than some threshold, let me know.

path_to_old_spectra='old/'
masks=['Q0110i1','Q0110i2']

text=input('Running this will overwrite redshift table files with those in {:s}. If you\'re sure you want to run it, type "yes": '.format(path_to_old_spectra))
if text != 'yes':
   import sys
   sys.exit()

from redshifting import redshift
from shutil import copy
from astropy.table import Table,Column,vstack,unique
from astropy.io import fits
import os
import numpy as np
import pdb
c=2.9979e5

# These are copied over from cubs_redshifting.py since I can't import that
def createSpec1Dfiles(mask,version='carpy'):
   print('Creating spec1D files')

   if version=='carpy':
      # try: #Stupid Python 2/3 compatibility
         # import cPickle as pickle
      # except:
         # import pickle
      # gdict=pickle.load(open('{}.p'.format(mask),'rb'))
      maskinfo=Table.read('{}/{}_maskinfo.fits'.format(mask,mask))
      path = '{}_spec1D'.format(mask)
      os.mkdir(path)
      # nRows=len(gdict['objects'])#-1 #last row is 'trace'
      nRows=len(maskinfo)#-1 #last row is 'trace'
      rows = Column(np.arange(nRows, dtype='int') + 1, name='row')
      ids = Column(np.chararray(nRows, itemsize=20), name='id')
      classes = Column(np.chararray(nRows, itemsize=6), name='class')
      redshifts = Column(np.zeros(nRows), name='redshift')
      qualities = Column(np.zeros(nRows, dtype='int') - 1, name='quality')
      comments = Column(np.chararray(nRows, itemsize=100), name='comment')
      extpos = Column(np.zeros(nRows,dtype='float'), name='extpos')
      extaper = Column(np.zeros(nRows,dtype='float'), name='extaper')
      extflag = Column(np.zeros(nRows,dtype=bool), name='extflag')
      boxes = Column(np.zeros(nRows,dtype=bool), name='alignbox')

      objects = Table([rows, ids, classes, redshifts, qualities, comments, extpos, extaper,extflag,boxes])
      objects['comment'] = 'none'
      objects['class'] = 'galaxy'
      objects['quality'] = -1
      # for i in range(nRows):
      for i,target in enumerate(maskinfo):
         # target=gdict['objects'][i]['setup']
         if os.path.isfile(target['spec1d']):
            spec1D=fits.getdata(target['spec1d'])
            header=fits.getheader(target['spec1d'])
            flux1Draw=spec1D[0,:] #No flux calibration
            error1Draw=spec1D[1,:]
            flux1D=spec1D[4,:] #Craptastic flux calibration
            error1D=spec1D[5,:]
            wave=(np.arange(header['NAXIS1'])+1-header['CRPIX1'])*header['CD1_1']+header['CRVAL1']
            if header['DC-FLAG']:
               wave=np.power(10,wave)
            #Convert wavelengths to vacuum from air
            #s = 1e4 / wave
            #n = (1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s**2)
            #     + 0.0001599740894897 / (38.92568793293 - s**2))
            #wave*=n
            spec = formatspectrum(wave, flux1D, error1D,raw=flux1Draw,rawerr=error1Draw)
            objID=target['object']
            objects[i]['id'] = objID
            if target['smf_type']=='HOLE': 
               objects[i]['class']='star'
               objects[i]['alignbox']=True
            objects[i]['extpos']=header['EXTRPOS']
            objects[i]['extaper']=header['APER']
         else:
            print('File is missing: {}'.format(target['spec1d']))
            objects[i]['id']=target['object']
            continue
         savename = getspec1Dname(mask, i+1, objID)
         fits.writeto(savename, spec)
      objects.write('{}/{}_objects.fits'.format(path, mask), overwrite=True)

def formatspectrum(wave, flux, error=0.0, mask=1, model=0.0, flat=0.0, arc=0.0,raw=0.0,rawerr=0.0):
      spec = np.zeros(len(flux),
                            dtype={'names':('wave', 'flux', 'error',
                                             'mask', 'model', 'flat', 'arc', 'raw', 'rawerr'), 
                                   'formats':(float, float, float, float, float,
                                              float, float,float,float)})
      spec['wave'] = wave
      spec['flux'] = flux
      spec['error'] = error
      spec['model'] = model
      spec['mask'] = mask
      spec['flat'] = flat
      spec['arc'] = arc
      spec['raw'] = raw
      spec['rawerr'] = rawerr
      
      
      spec['error'][~np.isfinite(spec['flux'])] = 0.0
      spec['flux'][~np.isfinite(spec['flux'])] = 0.0
      spec['flux'][~np.isfinite(spec['error'])] = 0.0
      spec['error'][~np.isfinite(spec['error'])] = 0.0
      
      index = np.where(spec['error'] == 0.0)
      spec['mask'][index] = 0.0
      
      return spec

def getspec1Dname(mask, row, id):
   
   return '{}_spec1D/{}_{}_{}.fits'.format(mask, mask, row, id)

def getspec2Dname(mask,id):

   return '{}/{}sum.fits'.format(mask,id)

def getredshift1Dname(mask, row, id):

   return '{}_spec1D/{}_{}_{}_redshift.fits'.format(mask, mask, row, id)

#Needed to modify this, since it assumes it's running inside the redshifting gui class
def redshiftObject(spec,obj,maskname):
      
      nGoodPix = np.sum(spec['mask'])
      if nGoodPix > 5:
         if obj['class'] == 'galaxy':
            redshifts = redshift.findz_galaxy(spec, zmin=-0.01, zmax=1.5, dz=0.0005)                     
            minIndex = np.argmin(redshifts['chi2_pdf'])
            z = redshifts[minIndex]['z']            
            redshifts_fine = redshift.findz_galaxy(spec, zmin=z-0.01, zmax=z+0.01, dz=0.0001)            
            redshifts = vstack((Table(redshifts), Table(redshifts_fine)))            
            redshifts = unique(redshifts, keys='z')            
            redshifts.sort('z')                        
            minIndex = np.argmin(redshifts['chi2_pdf'])
            z = redshifts[minIndex]['z']
            redshifts = np.array(redshifts)         
            eigenvalues, model, chi2pdf = redshift.fitatz_galaxy(spec, z)
            spec['model'] = model
            
         if obj['class'] == 'star':
            redshifts = redshift.findz_star(spec, zmin=-0.01, zmax=0.01, dz=0.0001)            
            minIndex = np.argmin(redshifts['chi2_pdf'])
            z = redshifts[minIndex]['z']
            eigenvalues, model, chi2pdf = redshift.fitatz_star(spec, z)
            spec['model'] = model  
            
            
         if obj['class'] == 'quasar':
            redshifts = redshift.findz_qso(spec, zmin=-0.01, zmax=4.0, dz=0.001)
            minIndex = np.argmin(redshifts['chi2_pdf'])
            z = redshifts[minIndex]['z']
            eigenvalues, model, chi2pdf = redshift.fitatz_qso(spec, z)
            spec['model'] = model
                                 
         obj['redshift'] = z
         redshift_savename = getredshift1Dname(maskname,obj['row'],obj['id'])
         Table(redshifts).write(redshift_savename,overwrite=True)


for mask in masks:
   print('\n{}'.format(mask))
   objfile=mask+'_objects.fits'
   objdir=mask+'_spec1D/'
   if not os.path.isdir(objdir):
      createSpec1Dfiles(mask)
   old_dir=path_to_old_spectra+objdir
   # old_dir=mask+'.old'
   copy(old_dir+objfile,objdir)
   objs=Table.read(objdir+objfile)
   for obj in objs:
      specname='{}_{:d}_{}.fits'.format(mask,obj['row'],obj['id'])
      oldspec=fits.getdata(old_dir+objdir+specname)
      newspec=fits.getdata(objdir+specname)
      newspec['mask']=oldspec['mask']
      ##UNCONNEMT THESE LINES TO UPDATE REDSHIFTS
      # if obj['quality']==2:
         # oldz=obj['redshift']
         # redshiftObject(newspec,obj,mask)
         # newz=obj['redshift']
         # dv=(newz-oldz)*c/(1+min([oldz,newz]))
         # print(obj['row'],obj['id'],oldz,newz,dv)
      fits.update(objdir+specname,newspec,1)
   objs.write(objdir+objfile,overwrite=True)

