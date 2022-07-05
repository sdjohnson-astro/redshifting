#!/usr/bin/env python
from PyQt5 import QtGui, QtCore  # (the example applies equally well to PySide)
import pyqtgraph as pg
import sys
import os
from astropy.io import fits
from astropy.table import Table, Column, vstack, unique
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from astropy.table import Table
import numpy as np
import glob
import argparse
import pyqtgraph.parametertree as pt
import redshift
import shutil




def formatspectrum(wave, flux, error=0.0, mask=1, model=0.0, flat=0.0, arc=0.0,raw=0.0,rawerr=0.0):
      spec = np.zeros(len(flux),
                            dtype={'names':('wave', 'flux', 'error',
                                             'mask', 'model', 'flat', 'arc', 'raw', 'rawerr'), 
                                   'formats':(float, float, float, float, float,
                                              float, float,float,float)})
      # spec['wave'] = 0.0
      # spec['flux'] = 0.0
      # spec['error'] = 0.0
      # spec['mask'] = 1
      # spec['model'] = 0.0
      # spec['flat'] = 0.0
      # spec['arc'] = 0.0  

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

#for CarPy output
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
   elif version=='cosmos':
      print('Creating spec1D files')
      spec1Darray = fits.getdata(mask + '_1spec.fits')
      header1D = fits.getheader(mask + '_1spec.fits')
      # Create wavelength array
      wave = header1D['CRVAL1'] + np.arange(header1D['NAXIS1'])*header1D['CDELT1']
      nRows = spec1Darray.shape[1]
      path = '{}_spec1D'.format(mask)
      os.mkdir(path)

      rows = Column(np.arange(nRows, dtype='int') + 1, name='row')
      ids = Column(np.chararray(nRows, itemsize=20), name='id')
      classes = Column(np.chararray(nRows, itemsize=6), name='class')
      redshifts = Column(np.zeros(nRows), name='redshift')
      qualities = Column(np.zeros(nRows, dtype='int') - 1, name='quality')
      comments = Column(np.chararray(nRows, itemsize=100), name='comment')
      #Things below here are just included for consistency with the carpy format...
      extpos = Column(np.zeros(nRows,dtype='float'), name='extpos')
      extaper = Column(np.zeros(nRows,dtype='float'), name='extaper')
      extflag = Column(np.zeros(nRows,dtype=bool), name='extflag')
      boxes = Column(np.zeros(nRows,dtype=bool), name='alignbox')
      # objects = Table([rows, ids, classes, redshifts, qualities, comments])
      objects = Table([rows, ids, classes, redshifts, qualities, comments, extpos, extaper,extflag])
      objects['comment'] = 'none'
      
      for i in range(nRows):
         flux1Draw = spec1Darray[0, i, :]
         error1Draw = spec1Darray[1, i, :]
         flat1D = spec1Darray[2, i, :]
         arc1D = spec1Darray[3, i, :]
         flux1D = spec1Darray[4, i, :]
         error1D = spec1Darray[5, i, :]
         spec = formatspectrum(wave, flux1D, error1D, 1.0, 0.0, flat1D, arc1D)
         apnum1D = header1D['APNUM{}'.format(i+1)]
         apnum1Darray = apnum1D.split(' ')
         id = apnum1Darray[1]
         savename = getspec1Dname(mask, i+1, id)
         fits.writeto(savename, spec)
         objects[i]['id'] = id
         print(i+1)         
      objects['class'] = 'galaxy'
      objects['quality'] = -1  
      objects.write('{}/{}_objects.fits'.format(path, mask), overwrite=True)         


class ldss3_redshiftgui:
   """Defining the GUI class for LDSS3 redshift assignment"""
   
   def __init__(self, mask, xsize=1000, ysize=1000, version='carpy'):
      
      self.mask = mask
      self.xsize = xsize
      self.ysize = ysize
      self.version = version
      #Read in the extraction profile for CarPy version
      if self.version=='carpy':
         try:
            # import pickle
            # gdict=pickle.load(open('{}.p'.format(mask),'rb'))
            # self.exttrace=gdict['trace'][0]
            self.exttrace=fits.getdata('{}/{}_trace.fits'.format(mask,mask))
         except:
            print('Couldn\'t read extraction trace file')
      
      # Read in the 1d and 2d files
      if self.version=='cosmos':
         self.spec1Darray = fits.getdata(mask + '_1spec.fits')
         self.header1D = fits.getheader(mask + '_1spec.fits')
         self.spec2Darray = fits.getdata(mask + '_big.fits')
         self.header2D = fits.getheader(mask + '_big.fits')
            
      path = '{}_spec1D'.format(mask)
      self.objects = Table.read('{}/{}_objects.fits'.format(path, mask))

      
      # Set the initial row number to zero
      self.row = 1
      self.nRows = len(self.objects)
      self.smoothing = 1
      
      self.z = 0.0
      self.redshifted = 0
      
      
      # masking flags
      self.mask_flag = 0
      self.mask_left =  np.array([0.0, 0.0])
      self.mask_right =  np.array([0.0, 0.0])
      
      
      
      # Get the GUI ready
      self.app = QtGui.QApplication([])
      #self.app = QtWidgets.QApplication(sys.argv)       # Always start by initializing Qt
      self.widget = QtGui.QWidget()       # Define a top-level widget
      
      # Set the widget size
      self.widget.resize(self.xsize, self.ysize)
      
      
      # Set the background plotting widget
      
      self.plot_redshift = pg.PlotWidget()
      self.plot_redshift.getAxis('bottom').setPen(pg.mkPen('w', width=2))
      self.plot_redshift.getAxis('top').setPen(pg.mkPen('w', width=2))
      self.plot_redshift.getAxis('left').setPen(pg.mkPen('w', width=2))
      self.plot_redshift.getAxis('right').setPen(pg.mkPen('w', width=2))
      self.plot_redshift.getAxis('bottom').setStyle(tickLength=-15)
      self.plot_redshift.getAxis('top').setStyle(tickLength=-15)
      self.plot_redshift.getAxis('left').setStyle(tickLength=-15)
      self.plot_redshift.getAxis('right').setStyle(tickLength=-15)
      self.plot_redshift.showAxis('right')
      self.plot_redshift.showAxis('top')
      self.plot_redshift.setLabel('bottom', 'redshift')
      self.plot_redshift.setLabel('left', 'chi2')
      
      
      
      
      # Set the 2D spectrum
      
      #self.plot_spec2D = pg.ImageView()
      #self.plot_spec2D.removeItem(self.plot_spec2D.getHistogramWidget())
      self.plot_spec2D_win = pg.GraphicsLayoutWidget()
      # self.plot_spec2D_view = self.plot_spec2D_win.addViewBox()
      self.plot_spec2D_plot = self.plot_spec2D_win.addPlot()      
      self.plot_spec2D = pg.ImageItem(border='w')
      self.plot_spec2D_plot.addItem(self.plot_spec2D)
      self.plot_spec2D_plot.setMouseEnabled(x=False, y=False)
      self.plot_spec2D_hist = pg.HistogramLUTWidget()
      self.plot_spec2D_hist.setImageItem(self.plot_spec2D)
      cm = self.plot_spec2D_hist.gradient.colorMap()
      cm.pos=np.array([1.,0.]) #is this really the easiest way to make white->black into black->white?
      self.plot_spec2D_hist.gradient.setColorMap(cm)
      
      # self.plot_spec2D.scene().sigMouseMoved.connect(self.mouseMoved_spec2D)
      
      # Set the 1D spectrum
      self.plot_spec1D = pg.PlotWidget()
      self.plot_spec1D.getAxis('bottom').setPen(pg.mkPen('w', width=2))
      self.plot_spec1D.getAxis('top').setPen(pg.mkPen('w', width=2))
      self.plot_spec1D.getAxis('left').setPen(pg.mkPen('w', width=2))
      self.plot_spec1D.getAxis('right').setPen(pg.mkPen('w', width=2))
      self.plot_spec1D.getAxis('bottom').setStyle(tickLength=-15)
      self.plot_spec1D.getAxis('top').setStyle(tickLength=-15)
      self.plot_spec1D.getAxis('left').setStyle(tickLength=-15)
      self.plot_spec1D.getAxis('right').setStyle(tickLength=-15)
      self.plot_spec1D.showAxis('right')
      self.plot_spec1D.showAxis('top')
      self.plot_spec1D.setLabel('bottom', 'Wavelength [&#8491;]')

      # self.plot_spec2D.getAxis('bottom').linkToView(self.plot_spec1D.getVieWBox())
      #self.plot_spec1D.setLabel('left', 'Flux')
      
      self.mouse_x_spec1D = 0.0
      self.mouse_y_spec1D = 0.0
      self.mouse_x_redshift = 0.0
      self.mouse_y_redshift = 0.0
      self.mouse_x_spec2D = 0.0
      self.mouse_y_spec2D = 0.0
      
      # Tie the image and spectrum
      # self.plot_spec2D.translate(4000,0)
      # self.plot_spec2D.scale(2,1)
      # self.plot_spec2D_view.linkView(self.plot_spec2D_view.XAxis,self.plot_spec1D.getViewBox())
      
      
      # Setup the layout
      self.layout = QtGui.QGridLayout()
      self.widget.setLayout(self.layout)
      
      
      
      # Set up right click menu
      self.featureListMenu = QtGui.QMenu("Galaxy lines")
      
      self.OVI1033 = QtGui.QAction("OVI 1033.82", self.featureListMenu)
      self.OVI1033.triggered.connect(self.setRedshiftOVI1033)
      self.featureListMenu.addAction(self.OVI1033)
      
      
      self.HI1215 = QtGui.QAction("HI Lya 1215.24", self.featureListMenu)
      self.HI1215.triggered.connect(self.setRedshiftHI1215)
      self.featureListMenu.addAction(self.HI1215)
      
      self.NV1240 = QtGui.QAction("NV 1240.81", self.featureListMenu)
      self.NV1240.triggered.connect(self.setRedshiftNV1240)
      self.featureListMenu.addAction(self.NV1240)
      
      self.CIV1549 = QtGui.QAction("CIV 1549.48", self.featureListMenu)
      self.CIV1549.triggered.connect(self.setRedshiftCIV1549)
      self.featureListMenu.addAction(self.CIV1549)
      
      self.OIII1665 = QtGui.QAction("OIII 1665.85", self.featureListMenu)
      self.OIII1665.triggered.connect(self.setRedshiftOIII1665)
      self.featureListMenu.addAction(self.OIII1665)
      
      self.CIII1908 = QtGui.QAction("CIII 1908.734", self.featureListMenu)
      self.CIII1908.triggered.connect(self.setRedshiftCIII1908)
      self.featureListMenu.addAction(self.CIII1908)
      
      self.MgII2799 = QtGui.QAction("MgII 2799.117", self.featureListMenu)
      self.MgII2799.triggered.connect(self.setRedshiftMgII2799)
      self.featureListMenu.addAction(self.MgII2799)
      
      self.OII3728 = QtGui.QAction("[OII] 3728.60", self.featureListMenu)
      self.OII3728.triggered.connect(self.setRedshiftOII3728)
      self.featureListMenu.addAction(self.OII3728)
      
      self.CaIIK3934 = QtGui.QAction("CaII K 3934.777", self.featureListMenu)
      self.CaIIK3934.triggered.connect(self.setRedshiftCaIIK3934)
      self.featureListMenu.addAction(self.CaIIK3934)
      
      self.CaIIH3969 = QtGui.QAction("CaII H 3969.588", self.featureListMenu)
      self.CaIIH3969.triggered.connect(self.setRedshiftCaIIH3969)
      self.featureListMenu.addAction(self.CaIIH3969)
      
      self.Hd4102 = QtGui.QAction("Hd 4102.89", self.featureListMenu)
      self.Hd4102.triggered.connect(self.setRedshiftHd4102)
      self.featureListMenu.addAction(self.Hd4102)
      
      self.Gband4305 = QtGui.QAction("G-band 4305.61", self.featureListMenu)
      self.Gband4305.triggered.connect(self.setRedshiftGband4305)
      self.featureListMenu.addAction(self.Gband4305)
      
      self.Hg4341 = QtGui.QAction("Hg 4341.68", self.featureListMenu)
      self.Hg4341.triggered.connect(self.setRedshiftHg4341)
      self.featureListMenu.addAction(self.Hg4341)
      
      self.OIII4364 = QtGui.QAction("[OIII] 4364.436", self.featureListMenu)
      self.OIII4364.triggered.connect(self.setRedshiftOIII4364)
      self.featureListMenu.addAction(self.OIII4364)
      
      self.Hb4862 = QtGui.QAction("Hb 4862", self.featureListMenu)
      self.Hb4862.triggered.connect(self.setRedshiftHb4862)
      self.featureListMenu.addAction(self.Hb4862)
      
      self.OIII4960 = QtGui.QAction("[OIII] 4960.295", self.featureListMenu)
      self.OIII4960.triggered.connect(self.setRedshiftOIII4960)
      self.featureListMenu.addAction(self.OIII4960)
      
      self.OIII5008 = QtGui.QAction("[OIII] 5008.240", self.featureListMenu)
      self.OIII5008.triggered.connect(self.setRedshiftOIII5008)
      self.featureListMenu.addAction(self.OIII5008)
      
      self.MgI5176 = QtGui.QAction("MgI 5176.7", self.featureListMenu)
      self.MgI5176.triggered.connect(self.setRedshiftMgI5176)
      self.featureListMenu.addAction(self.MgI5176)
      
      self.NaI5895 = QtGui.QAction("NaI 5895.6", self.featureListMenu)
      self.NaI5895.triggered.connect(self.setRedshiftNaI5895)
      self.featureListMenu.addAction(self.NaI5895)
      
      self.OI6302 = QtGui.QAction("[OI] 6302.046", self.featureListMenu)
      self.OI6302.triggered.connect(self.setRedshiftOI6302)
      self.featureListMenu.addAction(self.OI6302)
      
      self.OI6365 = QtGui.QAction("[OI] 6365.536", self.featureListMenu)
      self.OI6365.triggered.connect(self.setRedshiftOI6365)
      self.featureListMenu.addAction(self.OI6365)
      
      self.NII6549 = QtGui.QAction("[NII] 6549.86", self.featureListMenu)
      self.NII6549.triggered.connect(self.setRedshiftNII6549)
      self.featureListMenu.addAction(self.NII6549)
      
      self.Ha6564 = QtGui.QAction("Ha 6564.61", self.featureListMenu)
      self.Ha6564.triggered.connect(self.setRedshiftHa6564)
      self.featureListMenu.addAction(self.Ha6564)
      
      self.NII6585 = QtGui.QAction("[NII] 6585.27", self.featureListMenu)
      self.NII6585.triggered.connect(self.setRedshiftNII6585)
      self.featureListMenu.addAction(self.NII6585)
      
      self.SII6718 = QtGui.QAction("[SII] 6718.29", self.featureListMenu)
      self.SII6718.triggered.connect(self.setRedshiftSII6718)
      self.featureListMenu.addAction(self.SII6718)
      
      self.SII6732 = QtGui.QAction("[SII] 6732.67", self.featureListMenu)
      self.SII6732.triggered.connect(self.setRedshiftSII6732)
      self.featureListMenu.addAction(self.SII6732)
      
      self.CaII8500 = QtGui.QAction("CaII 8500.36", self.featureListMenu)
      self.CaII8500.triggered.connect(self.setRedshiftCaII8500)
      self.featureListMenu.addAction(self.CaII8500)
      
      self.CaII8544 = QtGui.QAction("CaII 8544.44", self.featureListMenu)
      self.CaII8544.triggered.connect(self.setRedshiftCaII8544)
      self.featureListMenu.addAction(self.CaII8544)
      
      self.CaII8664 = QtGui.QAction("CaII 8664.52", self.featureListMenu)
      self.CaII8664.triggered.connect(self.setRedshiftCaII8664)
      self.featureListMenu.addAction(self.CaII8664)
      
      
      
      self.plot_spec1D.getPlotItem().ctrlMenu = []
      self.plot_spec1D.getPlotItem().ctrlMenu = [self.featureListMenu]
      
      
      # Add listeners
      self.plot_spec1D.scene().sigMouseMoved.connect(self.mouseMoved_spec1D)
      self.plot_spec1D.keyPressEvent = self.keypress_spec1D
      
      self.plot_redshift.scene().sigMouseMoved.connect(self.mouseMoved_redshift)
      self.plot_redshift.keyPressEvent = self.keypress_redshift
      
      
      self.paramSpec = [
              dict(name='z=', type='str', value=self.z, dec=False, step=0.0001, limits=[None, None], readonly=True),
              dict(name='quality:', type='str', value='', readonly=True),
              dict(name='class:', type='str', value='', readonly=True), 
              dict(name='row:', type='str', value='', readonly=True),
              dict(name='id:', type='str', value='', readonly=True), 
              dict(name='Show lines', type='bool', value=True),
              dict(name='Show trace', type='bool', value=True),
              dict(name='Show raw', type='bool', value=False),
              dict(name='Bad Extraction:', type='bool', value=False)
              # dict(name='extraction center:', type='str', value='', readonly=True)
              # dict(name='extraction aper:', type='str', value='', readonly=True)
           ]
           
      self.param = pt.Parameter.create(name='Options', type='group', children=self.paramSpec)
      #Redraw when the boolean option buttons are pressed
      self.param.children()[5].sigValueChanged.connect(self.draw)
      self.param.children()[6].sigValueChanged.connect(self.draw)
      self.param.children()[7].sigValueChanged.connect(self.draw)
      self.param.children()[8].sigValueChanged.connect(self.setExtFlag)
      self.tree = pt.ParameterTree()
      self.tree.setParameters(self.param)

      
      self.features = Table.read(os.environ['REDSHIFTING'] + '/redshiftLines.dat', format='ascii')
      
      
      self.objectsTable = pg.TableWidget(editable=False, sortable=False)
      self.objectsTable.setFormat('%0.5f', 2)
      self.setTable()
      
      self.objectsTable.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
      self.objectsTable.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
      self.objectsTable.doubleClicked.connect(self.goToObject)
      
      # Add comment bar
      self.comment_text = QtGui.QLineEdit('comments here')
      self.comment_text.focusOutEvent = self.updateComment
      
      
      # Add plot_bg to the layout
      self.layout.addWidget(self.plot_redshift, 0, 0)
      self.layout.addWidget(self.plot_spec2D_win, 1, 0)
      self.layout.addWidget(self.plot_spec2D_hist, 1, 1)
      self.layout.addWidget(self.tree, 0, 1)
      self.layout.addWidget(self.objectsTable, 2, 1)
      self.layout.addWidget(self.plot_spec1D, 2, 0)
      self.layout.addWidget(self.comment_text, 4, 0)
      
      
      self.layout.setColumnStretch(0, 4)
      self.layout.setColumnStretch(1, 1)
      self.layout.setColumnStretch(1, 1)
      
      
      
      #self.layout.setColumnStretch(0, 2)
      #self.layout.setColumnStretch(0, 2)
      #self.layout.setColumnStretch(3, 1)
      #self.layout.setColumnStretch(3, 0)
      self.layout.setRowStretch(0, 2)
      self.layout.setRowStretch(1, 1)
      self.layout.setRowStretch(2, 3)
      
      
      self.setSpec()

      #Set 2D X-axis values to be equal to the 1D--this works so long as all spectra are the same size
      #since it is just based on the first one read in
      self.plot_spec2D.translate(min(self.wave),0)
      self.plot_spec2D.scale((max(self.wave)-min(self.wave))/len(self.wave),1)
      self.plot_spec2D_plot.setXLink(self.plot_spec1D)

      
      self.draw()
      
      self.widget.show()
      self.app.exec_()
      
   
   def updateComment(self, event):
      
      self.objects[self.row-1]['comment'] = self.comment_text.text()
         
   def setTable(self):
      
      self.objectsTable.setData(np.array(self.objects['id', 'class', 'redshift', 'quality', 'extflag','row','comment']))
   
   def goToObject(self):
      
      #print('Going to object...')
      self.row = self.objectsTable.selectedItems()[0].row()+1
      self.setSpec()
      # self.draw()
   
   def setRedshiftOVI1033(self):
         wave0 = 1033.82
         self.z = self.mouse_x_spec1D/wave0 - 1
         #self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
            
   
   def setRedshiftHI1215(self):
         wave0 = 1215.24
         self.z = self.mouse_x_spec1D/wave0 - 1
         #self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftNV1240(self):
         wave0 = 1240.81
         self.z = self.mouse_x_spec1D/wave0 - 1
         #self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftCIV1549(self):
         wave0 = 1549.48
         self.z = self.mouse_x_spec1D/wave0 - 1
         #self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftOIII1665(self):
         wave0 = 1665.85
         self.z = self.mouse_x_spec1D/wave0 - 1
         #self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
            
   def setRedshiftCIII1908(self):
         wave0 = 1908.734
         self.z = self.mouse_x_spec1D/wave0 - 1
         #self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftMgII2799(self):
         wave0 = 2799.117
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftOII3728(self):
         wave0 = 3728.60
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftCaIIK3934(self):
         wave0 = 3934.777
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftCaIIH3969(self):
         wave0 = 3969.588
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftHd4102(self):
         wave0 = 4102.89
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftGband4305(self):
         wave0 = 4305.61
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftHg4341(self):
         wave0 = 4341.68
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftOIII4364(self):
         wave0 = 4364.436
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftHb4862(self):
         wave0 = 4862.68
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftOIII4960(self):
         wave0 = 4960.295
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftOIII5008(self):
         wave0 = 5008.240
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftMgI5176(self):
         wave0 = 5176.7
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftNaI5895(self):
         wave0 = 5895.6
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftOI6302(self):
         wave0 = 6302.046
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftOI6365(self):
         wave0 = 6365.536
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftNII6549(self):
         wave0 = 6549.86
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftHa6564(self):
         wave0 = 6564.61
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.fitObjectAtRedshift()
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftNII6585(self):
         wave0 = 6585.27
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftSII6718(self):
         wave0 = 6718.29
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftSII6732(self):
         wave0 = 6732.67
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftCaII8500(self):
         wave0 = 8500.36
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftCaII8544(self):
         wave0 = 8544.44
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   def setRedshiftCaII8664(self):
         wave0 = 8664.52
         self.z = self.mouse_x_spec1D/wave0 - 1
         self.fitObjectAtRedshift()
         self.objects[self.row-1]['redshift'] = self.z
         self.param['z='] =  '{:0.5f}'.format(self.z)
         self.draw()
                  
   
   def setClass(self, classification):
         self.objects[self.row-1]['class'] = classification
         self.param['class:'] =  classification
         
         self.draw()
               

   def setQuality(self, quality):
      
      self.objects[self.row-1]['quality'] = quality
      self.param['quality:'] =  quality
      
      self.draw()
      
   def setExtFlag(self):
    
      self.objects[self.row-1]['extflag']=self.param['Bad Extraction:']
      self.save()
      
   def keypress_redshift(self, event):
      
      if self.plot_redshift_current: 
         
         #print('')
         #print(event.text())
         #print(event.key())
         #print('{:0.4f}, {:0.2f}'.format(self.mouse_x_redshift, self.mouse_y_redshift))
   
         if (event.text() == 'n'):
            
            self.advance(1)
         
         if (event.text() == 'N'):
            
            self.advance(10)
            
         if (event.text() == 'b'):
            
            self.advance(-1)
            
         if (event.text() == 'B'):
            
            self.advance(-10)
         
         if event.text() == '[':
            self.panx_redshift(-1.0/3.0)
            
         if event.text() == ']':
            self.panx_redshift(1.0/3.0)
            
         if event.text() == '{':
            self.panx_redshift(-1.0)
            
         if event.text() == '}':
            self.panx_redshift(1.0)
            
         if event.text() == 'x':
            self.zoomxy_redshift(1.0/1.5, 1.0)
         
         if event.text() == 'X':
            self.zoomxy_redshift(1.5, 1.0)
            
         if event.text() == 'y':
            self.zoomxy_redshift(1.0, 1.0/1.5)
         
         if event.text() == 'Y':
            self.zoomxy_redshift(1.0, 1.5)
            
         if (event.text() == 'w') | (event.text() == 'W'):
            self.plot_redshift.autoRange()
   
   
         # Set redshift
         if event.text() == 'z':
            self.z = self.mouse_x_redshift
            self.objects[self.row-1]['redshift'] = self.z
            self.param['z='] =  '{:0.5f}'.format(self.z)
            self.fitObjectAtRedshift()
            
            self.draw()
            
                        
         if (event.text() == 'h') | (event.text() == 'H'):
            
            self.setClass('star')
            
         if (event.text() == 'g') | (event.text() == 'G'):
            
            self.setClass('galaxy')
            
         if (event.text() == 'j') | (event.text() == 'J'):
            
            self.setClass('quasar')
            
         if (event.text() == 'k') | (event.text() == 'K'):
            
            self.setClass('hizgal')
            
            
         # if event.text() == 'R':
            
            # self.redshiftAll()
            
         
         if event.text() == 'r':
            
            self.redshiftObject()
            
         if event.text() == 'l':
            
            self.redshiftObjectLocal()
            
         if (event.text() == 'a') | (event.text() == 'A'):
            
            self.setQuality(2)
            
         if (event.text() == 's') | (event.text() == 'S'):
            
            self.setQuality(1)
         
         if (event.text() == 'd') | (event.text() == 'D'):
            
            self.setQuality(0)
            
         if (event.text() == 'f') | (event.text() == 'F'):
            
            self.setQuality(-1)
            
         if event.text() == ';':
            
            self.incrementRedshift(-0.0001)
            
         if event.text() == "'":
            
            self.incrementRedshift(+0.0001)
            
            
         if event.text() == ':':
            
            self.incrementRedshift(-0.001)
            
         if event.text() == '"':
            
            self.incrementRedshift(+0.001)
            
         if event.text() == ',':
            
            self.incrementRedshift(-0.01)
            
         if event.text() == '.':
            
            self.incrementRedshift(+0.01)
            
         if event.text() == '<':
            
            self.incrementRedshift(-0.1)
            
         if event.text() == '>':
            
            self.incrementRedshift(+0.1)
            
         # left arrow
         if event.key() == 16777234:
            
            self.incrementRedshift(-1)
         
         # Right arrow   
         if event.key() == 16777236:
            
            self.incrementRedshift(+1)            
            
   def zoomxy_redshift(self, scalex, scaley):
      """Zoom in or out in wavelength (x) and/or flux (y)"""
   
      xRange = self.plot_redshift.getViewBox().state['viewRange'][0]
      x0 = xRange[0]
      x1 = xRange[1]
      xRange = (x1 - x0)*scalex
      x0_new = self.mouse_x_redshift - xRange/2.0
      x1_new = self.mouse_x_redshift + xRange/2.0
      
      self.plot_redshift.setXRange(x0_new, x1_new, padding=0)
      
      yRange = self.plot_redshift.getViewBox().state['viewRange'][1]
      y0 = yRange[0]
      y1 = yRange[1]
      yRange = (y1 - y0)*scaley
      y0_new = self.mouse_y_redshift - yRange/2.0
      y1_new = self.mouse_y_redshift + yRange/2.0
      
      self.plot_redshift.setYRange(y0_new, y1_new, padding=0)

   def zoom_default(self):
      q1,q2=np.percentile(self.flux1D[self.spec['mask'].astype('bool')],[1,99])
      q1 = np.min([0.0, q1])
      self.plot_spec1D.setYRange(-0.05*q2,2*q2,padding=0)
      self.plot_spec1D.setXRange(np.min(self.wave), np.max(self.wave), padding=0)

   
   # def zoom_redshift(self, key, ):
   #    """Zoom in or out in wavelength (x) and/or flux (y)"""
   
   #    xRange = self.plot_redshift.getViewBox().state['viewRange'][0]
   #    x0 = xRange[0]
   #    x1 = xRange[1]
   #    xRange = (x1 - x0)*scalex
   #    x0_new = self.mouse_x_redshift - xRange/2.0
   #    x1_new = self.mouse_x_redshift + xRange/2.0
      
   #    self.plot_redshift.setXRange(x0_new, x1_new, padding=0)
      
   #    yRange = self.plot_redshift.getViewBox().state['viewRange'][1]
   #    y0 = yRange[0]
   #    y1 = yRange[1]
   #    yRange = (y1 - y0)*scaley
   #    y0_new = self.mouse_y_redshift - yRange/2.0
   #    y1_new = self.mouse_y_redshift + yRange/2.0
      
   #    self.plot_redshift.setYRange(y0_new, y1_new, padding=0)
      
   def panx_redshift(self, scalex):
      """Pan in the wavelength direction"""
      
      xRange = self.plot_redshift.getViewBox().state['viewRange'][0]
      x0 = xRange[0]
      x1 = xRange[1]
      shift = scalex*(x1 - x0)
      x0_new = x0 + shift
      x1_new = x1 + shift

      self.plot_redshift.setXRange(x0_new, x1_new, padding=0)
            
            
   
   def keypress_spec1D(self, event):
      
      
      if self.plot_spec1D_current:
         
         #print(event.text())
         #print(event.key())
         #print('')
         
         if (event.text() == 'n'):
            
            self.advance(1)
         
         if (event.text() == 'N'):
            
            self.advance(10)
            
         if (event.text() == 'b'):
            
            self.advance(-1)
            
         if (event.text() == 'B'):
            
            self.advance(-10)
         
         if event.text() == '[':
            self.panx_spec(-1.0/3.0)
            
         if event.text() == ']':
            self.panx_spec(1.0/3.0)
            
         if event.text() == '{':
            self.panx_spec(-1.0)
            
         if event.text() == '}':
            self.panx_spec(1.0)
            
         if event.text() == 'x':
            self.zoomxy_spec(1.0/1.5, 1.0)
         
         if event.text() == 'X':
            self.zoomxy_spec(1.5, 1.0)
            
         if event.text() == 'y':
            self.zoomxy_spec(1.0, 1.0/1.5)
         
         if event.text() == 'Y':
            self.zoomxy_spec(1.0, 1.5)

         if event.text() in ['t','T','e','E']:
            #The not is because empty strings (pressing shift) get included too!
            self.trimxy_spec(event.text())
            
         if (event.text() == 'w') | (event.text() == 'W'):
            self.zoom_default()
            # self.plot_spec1D.autoRange()
            # self.updateXrange_1D()            
            
         if (event.text() == '=') | (event.text() == '+'):
            self.smoothing = self.smoothing + 2
            
            if self.smoothing == 3:
               self.smoothing = 5
                  
            self.smoothSpec()
            
         if (event.text() == '-') | (event.text() == '_'):
            self.smoothing = self.smoothing - 2
            
            if self.smoothing < 5:
               self.smoothing = 1
         
            self.smoothSpec()
            
         if event.text() == 'm':
            
            self.changeMask(0)
            
         if event.text() == 'M':
            
            # self.changeMask(0)
            self.autoMask()
            
         if event.text() == 'u':
            
            self.changeMask(1)
            
         if event.text() == 'U':
            
            self.changeMask(1)
            
         if (event.text() == 'h') | (event.text() == 'H'):
            
            self.setClass('star')
            
         if (event.text() == 'g') | (event.text() == 'G'):
            
            self.setClass('galaxy')
            
         if (event.text() == 'j') | (event.text() == 'J'):
            
            self.setClass('quasar')
            
         if (event.text() == 'k') | (event.text() == 'K'):
            
            self.setClass('hizgal')
            
            
         if event.text() == '/':
            
             self.redshiftAll()
            
         
         if event.text() == 'r':
            
            self.redshiftObject()
            
         if event.text() == 'l':
            
            self.redshiftObjectLocal()
            
         if (event.text() == 'a') | (event.text() == 'A'):
            
            self.setQuality(2)
            
         if (event.text() == 's') | (event.text() == 'S'):
            
            self.setQuality(1)
         
         if (event.text() == 'd') | (event.text() == 'D'):
            
            self.setQuality(0)
            
         if (event.text() == 'f') | (event.text() == 'F'):
            
            self.setQuality(-1)

         # if (event.text() == 'c') & (self.version=='carpy'):
            
            # self.setExtFlag()
            
         if event.text() == ';':
            
            self.incrementRedshift(-0.0001)
            
         if event.text() == "'":
            
            self.incrementRedshift(+0.0001)
            
            
         if event.text() == ':':
            
            self.incrementRedshift(-0.001)
            
         if event.text() == '"':
            
            self.incrementRedshift(+0.001)
            
         if event.text() == ',':
            
            self.incrementRedshift(-0.01)
            
         if event.text() == '.':
            
            self.incrementRedshift(+0.01)
            
         if event.text() == '<':
            
            self.incrementRedshift(-0.1)
            
         if event.text() == '>':
            
            self.incrementRedshift(+0.1)
            
         # left arrow
         if event.key() == 16777234:
            
            self.incrementRedshift(-1)
         
         # Right arrow   
         if event.key() == 16777236:
            
            self.incrementRedshift(+1)
            
            
         if event.text() == '1':
            
            self.setRedshiftHa6564()
            
         
         if event.text() == '2':
            
            self.setRedshiftOIII5008()
            
         
         if event.text() == '3':
            
            self.setRedshiftOIII4960()
            
         if event.text() == '4':
            
            self.setRedshiftHb4862()
            
         if event.text() == '5':
            
            self.setRedshiftHg4341()
            
         if event.text() == '6':
            
            self.setRedshiftHd4102()
         
         if event.text() == '7':
            
            self.setRedshiftOII3728()
            
         if event.text() == '8':
            
            self.setRedshiftCIII1908()
            
         if event.text() == '9':
            
            self.setRedshiftCIV1549()
            
         if event.text() == '0':
            
            self.setRedshiftHI1215()
            
         if event.text() == '!':
            
            self.setRedshiftNaI5895()
            
         if event.text() == '@':
            
            self.setRedshiftMgI5176()
            
         if event.text() == '#':
            
            self.setRedshiftGband4305()
            
         if event.text() == '$':
            
            self.setRedshiftCaIIH3969()
            
         if event.text() == '%':
            
            self.setRedshiftCaIIK3934()
            
         if event.text() == '^':
            
            self.setRedshiftMgII2799()
            
         if event.text() == 'z':
            
            self.fitObjectAtRedshift()
            
            
            
         # if return or enter is pressed then save.
         if (event.key() == 16777220) | (event.key() == 16777221) | (event.key() == 96) | (event.key() == 126):
            
            self.save()
   
   
   def incrementRedshift(self, dz):
      
      self.z = self.z + dz
      self.objects[self.row-1]['redshift'] = self.z
      self.param['z='] =  '{:0.5f}'.format(self.z)
      self.draw()
         
   def redshiftAll(self):
      
      # Start at row 1
      for object in self.objects:
         # if it already has a redshift:
            # continue
         self.row = object['row']
         self.setSpec()
         if self.redshifted:
            print('{}/{}   {}  already redshifted. Skipping'.format(self.row, self.nRows, object['class']))
         else:
            self.autoMask()            
            self.redshiftObject()
            print('{}/{}   {}   z_best={:0.4f}'.format(self.row, self.nRows, object['class'], object['redshift']))
         
      
   
   def fitObjectAtRedshift(self):
      
      spec = self.spec
      if self.objects[self.row-1]['class'] == 'galaxy':
         
         z = self.z
         eigenvalues, model, chi2pdf = redshift.fitatz_galaxy(spec, z)
         spec['model'] = model
         
         
      if self.objects[self.row-1]['class'] == 'star':
         
         z = self.z
         eigenvalues, model, chi2pdf = redshift.fitatz_star(spec, z)
         spec['model'] = model  
         
         
      if self.objects[self.row-1]['class'] == 'quasar':
         
         z = self.z
         eigenvalues, model, chi2pdf = redshift.fitatz_qso(spec, z)
         spec['model'] = model
         
      if self.objects[self.row-1]['class'] == 'hizgal':
         
         z = self.z
         eigenvalues, model, chi2pdf = redshift.fitatz_latis(spec, z)
         spec['model'] = model
         
         
      print('Redshift assigned by hand {}   {}   z={:0.4f} and saved'.format(self.objects[self.row-1]['row'],self.objects[self.row-1]['id'], z))

      
      self.objects[self.row-1]['redshift'] = z
      #self.redshifted = 1
      #self.redshifts = Table(redshifts)
      self.save()
      self.draw()      
   
   
   
   
   def redshiftObjectLocal(self):
      
      spec = self.spec
      
      
      z = self.z
      
      if self.objects[self.row-1]['class'] == 'galaxy':
         
         redshifts = redshift.findz_galaxy(spec, zmin=z-0.01, zmax=z+0.01, dz=0.0001)
                  
         minIndex = np.argmin(redshifts['chi2_pdf'])
         z = redshifts[minIndex]['z']
         
         eigenvalues, model, chi2pdf = redshift.fitatz_galaxy(spec, z)
         spec['model'] = model
         
         
      if self.objects[self.row-1]['class'] == 'star':
         
         redshifts = redshift.findz_star(spec, zmin=self.z-0.001, zmax=self.z+0.001, dz=0.0001)
         
         minIndex = np.argmin(redshifts['chi2_pdf'])
         z = redshifts[minIndex]['z']
         eigenvalues, model, chi2pdf = redshift.fitatz_star(spec, z)
         spec['model'] = model  
         
         
      if self.objects[self.row-1]['class'] == 'quasar':
         
         redshifts = redshift.findz_qso(spec, zmin=self.z-0.01, zmax=self.z+0.01, dz=0.001)
         
         minIndex = np.argmin(redshifts['chi2_pdf'])
         z = redshifts[minIndex]['z']
         eigenvalues, model, chi2pdf = redshift.fitatz_qso(spec, z)
         spec['model'] = model
         
      if self.objects[self.row-1]['class'] == 'hizgal':
         
         redshifts = redshift.findz_latis(spec, zmin=z-0.01, zmax=z+0.01, dz=0.0001)
                  
         minIndex = np.argmin(redshifts['chi2_pdf'])
         z = redshifts[minIndex]['z']
         
         eigenvalues, model, chi2pdf = redshift.fitatz_latis(spec, z)
         spec['model'] = model
      
         
      self.objects[self.row-1]['redshift'] = z
      if self.redshifted == 0:  
         
         self.redshifted = 1
         self.redshifts = Table(redshifts)
         
      else:
         
         redshifts = vstack((Table(self.redshifts), Table(redshifts)))
         redshifts = unique(redshifts, keys='z')
         redshifts.sort('z')
         self.redshifts = redshifts
         
      print('Redshifting Locally {}   {}   z={:0.4f} and saved'.format(self.objects[self.row-1]['row'],self.objects[self.row-1]['id'], z))

      self.z = z
      self.param['z='] =  '{:0.5f}'.format(self.z)       
      self.save()
      self.draw()
      
         
   def redshiftObject(self):
      
      spec = self.spec
      
      nGoodPix = np.sum(spec['mask'])
      if nGoodPix > 5:
      
      
         if self.objects[self.row-1]['class'] == 'galaxy':
            
            redshifts = redshift.findz_galaxy(spec, zmin=-0.01, zmax=1.5, dz=0.0005)
                     
            minIndex = np.argmin(redshifts['chi2_pdf'])
            z = redshifts[minIndex]['z']
            
            redshifts_fine = redshift.findz_galaxy(spec, zmin=z-0.003, zmax=z+0.003, dz=0.0001)
            
            redshifts = vstack((Table(redshifts), Table(redshifts_fine)))
            
            redshifts = unique(redshifts, keys='z')
            
            redshifts.sort('z')
            
            
            minIndex = np.argmin(redshifts['chi2_pdf'])
            z = redshifts[minIndex]['z']
            
            redshifts = np.array(redshifts)
            
            
            
            eigenvalues, model, chi2pdf = redshift.fitatz_galaxy(spec, z)
            spec['model'] = model
            
            
         if self.objects[self.row-1]['class'] == 'star':
            
            redshifts = redshift.findz_star(spec, zmin=-0.01, zmax=0.01, dz=0.0001)
            
            minIndex = np.argmin(redshifts['chi2_pdf'])
            z = redshifts[minIndex]['z']
            eigenvalues, model, chi2pdf = redshift.fitatz_star(spec, z)
            spec['model'] = model  
            
            
         if self.objects[self.row-1]['class'] == 'quasar':
            
            redshifts = redshift.findz_qso(spec, zmin=-0.01, zmax=4.0, dz=0.001)
            
            minIndex = np.argmin(redshifts['chi2_pdf'])
            z = redshifts[minIndex]['z']
            eigenvalues, model, chi2pdf = redshift.fitatz_qso(spec, z)
            spec['model'] = model
            
         
         if self.objects[self.row-1]['class'] == 'hizgal':
            
            redshifts = redshift.findz_latis(spec)
            
            minIndex = np.argmin(redshifts['chi2_pdf'])
            print(Table(redshifts))
            z = redshifts[minIndex]['z']
            eigenvalues, model, chi2pdf = redshift.fitatz_latis(spec, z)
            spec['model'] = model
            
            
         
         print('Redshifting {}   {}   z={:0.4f} and saved'.format(self.objects[self.row-1]['row'],
                                                        self.objects[self.row-1]['id'], z))
         
         self.objects[self.row-1]['redshift'] = z
         self.z = z
         self.redshifted = 1
         self.redshifts = Table(redshifts)
         self.param['z='] =  '{:0.5f}'.format(self.z)      
         self.save()
         self.draw()
         
      else:
         
         self.z = 0.0
         self.redshift = 0
         
      
   
            
   def save(self):
      
      self.setTable()
# 
      path = '{}_spec1D'.format(self.mask)
      self.objects.write('{}/{}_objects.fits'.format(path, self.mask), overwrite=True)

      savename = getspec1Dname(self.mask, self.objects[self.row-1]['row'],
                               self.objects[self.row-1]['id'])
      if not os.path.isfile(savename):
         return
      fits.writeto(savename, self.spec, overwrite=True)
      
      # If we have a redshift array, store it
      if self.redshifted == 1:
         savename = getredshift1Dname(self.mask, self.objects[self.row-1]['row'],
                                  self.objects[self.row-1]['id'])
         self.redshifts.write(savename, overwrite=True)
         
      print('Saved')   
      # if os.path.isfile(path):
      #    savename = getspec1Dname(self.mask, self.objects[self.row-1]['row'],
      #                          self.objects[self.row-1]['id'])
      #    fits.writeto(savename, self.spec, overwrite=True)
      
      # self.objects.write('{}/{}_objects.fits'.format(path, self.mask), overwrite=True)
      
      
      # # If we have a redshift array, store it
      # if self.redshifted == 1:
      #    savename = getredshift1Dname(self.mask, self.objects[self.row-1]['row'],
      #                             self.objects[self.row-1]['id'])
      #    self.redshifts.write(savename, overwrite=True)
         
      # print('Saved')   

   def autoMask(self):
      """Automatically mask things below and above some range, and the A band"""
      if (self.spec['wave'][0]<4000): #LDSS3
         low_cut=5000
         high_cut=9800
      else: #IMACS
         low_cut=5200
         high_cut=9275
      aband=[7588,7684]

      index = np.where((self.wave < low_cut) | (self.wave > high_cut)
         | ((self.wave > aband[0]) & (self.wave < aband[1])) | (self.error1D == 0))
      self.spec['mask'][index] = 0
      self.draw()

   def changeMask(self, newvalue):
      """Change the mask"""
      
      if self.mask_flag == 0:
         
         self.mask_left = np.array([self.mouse_x_spec1D, self.mouse_y_spec1D])
         
         self.mask_flag = 1
      
      else:
         self.mask_right = np.array([self.mouse_x_spec1D, self.mouse_y_spec1D])
         
         wave0 = np.min([self.mask_left[0], self.mask_right[0]])
         wave1 = np.max([self.mask_left[0], self.mask_right[0]])
         
         index = np.where((self.wave > wave0) & (self.wave < wave1))
         self.spec['mask'][index] = newvalue
         
         
         self.mask_flag = 0
         self.draw()
         
            
   def smoothSpec(self):
      """Smooth the spectrum using Savitzky-Golay filter."""
      
      if self.smoothing > 1:
         self.flux1D = savgol_filter(self.spec['flux'], self.smoothing, 2)
         self.error1D = savgol_filter(self.spec['error'], self.smoothing, 2)/np.sqrt(self.smoothing)
         self.model1D = savgol_filter(self.spec['model'], self.smoothing, 2)
      if self.smoothing == 1:
         self.flux1D = self.spec['flux']
         self.error1D = self.spec['error']
         self.model1D = self.spec['model']
      
      self.draw()

   def trimxy_spec(self, key):
      """Trim plotting region in wavelength (x) and/or flux (y)"""
      if key in 'eE':
         xRange = self.plot_spec1D.getViewBox().state['viewRange'][0]
         x0 = xRange[0]
         x1 = xRange[1]
         xnew = self.mouse_x_spec1D
         if key=='E': x1=xnew
         else: x0=xnew
         self.plot_spec1D.setXRange(x0, x1, padding=0)
      elif key in 'tT':
         yRange = self.plot_spec1D.getViewBox().state['viewRange'][1]
         y0 = yRange[0]
         y1 = yRange[1]
         ynew = self.mouse_y_spec1D
         if key=='t': y1=ynew
         else: y0=ynew
         self.plot_spec1D.setYRange(y0, y1, padding=0)
      self.updateXrange_1D()

            
   def zoomxy_spec(self, scalex, scaley):
      """Zoom in or out in wavelength (x) and/or flux (y)"""
   
      xRange = self.plot_spec1D.getViewBox().state['viewRange'][0]
      x0 = xRange[0]
      x1 = xRange[1]
      xRange = (x1 - x0)*scalex
      x0_new = self.mouse_x_spec1D - xRange/2.0
      x1_new = self.mouse_x_spec1D + xRange/2.0
      
      self.plot_spec1D.setXRange(x0_new, x1_new, padding=0)
      
      yRange = self.plot_spec1D.getViewBox().state['viewRange'][1]
      y0 = yRange[0]
      y1 = yRange[1]
      yRange = (y1 - y0)*scaley
      y0_new = self.mouse_y_spec1D - yRange/2.0
      y1_new = self.mouse_y_spec1D + yRange/2.0
      
      self.plot_spec1D.setYRange(y0_new, y1_new, padding=0)
      
      self.updateXrange_1D()

      # self.specCursor.setPos((x1_new-x0_new)/2,(y1_new-y0_new)/2)
      
   def panx_spec(self, scalex):
      """Pan in the wavelength direction"""
      
      xRange = self.plot_spec1D.getViewBox().state['viewRange'][0]
      x0 = xRange[0]
      x1 = xRange[1]
      shift = scalex*(x1 - x0)
      x0_new = x0 + shift
      x1_new = x1 + shift

      self.plot_spec1D.setXRange(x0_new, x1_new, padding=0)
      self.updateXrange_1D()
   
   
   
   def mouseMoved_spec1D(self, pos):
       """Keep track of where the mouse and update the xrange on the 2D plot to match the 1D"""
       self.plot_spec1D_current = True
       self.plot_spec2D_current = False
       self.plot_redshift_current = False
       
       self.mouse_x_spec1D = self.plot_spec1D.mapToView(pos).x()
       self.mouse_y_spec1D = self.plot_spec1D.mapToView(pos).y()
       
       self.updateXrange_1D()
       
       self.setTitle_1D()
       
    
       
   #This is now obsolete since the axes are linked
   def updateXrange_1D(self):
     pass 
      # xRange = self.plot_spec1D.getViewBox().state['viewRange'][0]
      # x0 = xRange[0]
      # x1 = xRange[1]
      
      # # Interpolate wavelength to index very inefficient way to do this but I am lazy
      # indexes = np.arange(len(self.wave))
      # indexes_interp = interp1d(self.wave, indexes, bounds_error=False, fill_value='extrapolate')
      # index0 = indexes_interp(x0)
      # index1 = indexes_interp(x1)
      
      # #index = np.where((self.wave > x0) & (self.wave < x1))[0]
      # #if len(index) > 0:
      # #index0 = np.min(index)
      # #index1 = np.max(index)
      # self.plot_spec2D_view.setXRange(index0, index1, padding=0.035)      
   
   
   def advance(self, delt):
      
      self.save()
      self.row = self.row + delt
      if self.row < 1:
         self.row = self.nRows
      if self.row > self.nRows:
         self.row = 1
      self.setSpec()
               
      #print('{}/{}'.format(self.row, self.nRows))
      #self.draw()
      
   def setSpec(self, autoRange=True):
      """Set the spectrum to current row"""
      self.z = self.objects[self.row-1]['redshift']
      self.comment_text.setText(self.objects[self.row-1]['comment'])

      if self.version=='carpy':
         self.id=self.objects[self.row-1]['id']
         if not os.path.isfile(getspec1Dname(self.mask, self.row, self.id)):
            print('Files for {} do not exist. Moving on.'.format(self.id))
            self.redshifted=0
            self.advance(1)
            return
         self.flux2D=fits.getdata(getspec2Dname(self.mask,self.id)).transpose()

      elif self.version=='cosmos':
         # Get the apnum header parameter
         self.apnum1D = self.header1D['APNUM{}'.format(self.row)]
         self.apnum2D = self.header2D['APNUM{}'.format(self.row)]
         self.apnum1Darray = self.apnum1D.split(' ')
         self.apnum2Darray = self.apnum2D.split(' ')
         self.id = self.apnum1Darray[1]
         self.y0=int(self.header2D['CSECT{}A'.format(self.row)])
         self.y1=int(self.header2D['CSECT{}B'.format(self.row)])
         self.y0 = int(float(self.apnum2Darray[2]))
         self.y1 = int(float(self.apnum2Darray[3]))
         self.flux2D = self.spec2Darray[self.y0:self.y1, :].transpose()
      
      self.spec = fits.getdata(getspec1Dname(self.mask, self.row, self.id))
      
      self.wave = self.spec['wave']
      self.flux1D = self.spec['flux']
      self.flux1Draw = self.spec['raw']
      self.error1D = self.spec['error']
      self.error1Draw=self.spec['rawerr']
      self.model1D = self.spec['model']
      self.flat1D = self.spec['flat']
      self.arc1D = self.spec['arc']

      self.extpos = self.objects[self.row-1]['extpos']
      self.extaper = self.objects[self.row-1]['extaper']
      self.smoothSpec()
      
      # Check for redshift filename and read in if present.
      redshiftFilename = getredshift1Dname(self.mask, self.objects[self.row-1]['row'],
                                           self.objects[self.row-1]['id'])
      if os.path.isfile(redshiftFilename):
         self.redshifted = 1
         self.redshifts = Table.read(redshiftFilename)
         #print('Already redshifted')
         
      else:
         self.redshifted = 0
         self.redshifts = None
      
      self.draw()
      self.plot_redshift.autoRange()
      if autoRange:
         # self.plot_spec1D.autoRange()
         self.zoom_default()
      self.plot_spec2D_hist.setHistogramRange(*np.percentile(self.flux2D,[0.1,99.9]))
      self.plot_spec2D_hist.setLevels(*np.percentile(self.flux2D,[0.5,99.5]))
      
      self.param['row:'] =  '{}/{}'.format(self.row, self.nRows)
      self.param['id:'] =  '{}'.format(self.objects[self.row-1]['id'])
      self.param['class:'] =  '{}'.format(self.objects[self.row-1]['class'])
      self.param['z='] =  '{:0.5f}'.format(self.objects[self.row-1]['redshift'])
      self.param['quality:'] =  '{}'.format(self.objects[self.row-1]['quality'])
      self.param['Bad Extraction:'] =  bool(self.objects[self.row-1]['extflag'])

      
   def mouseMoved_redshift(self, pos):
       """Keep track of where the mouse and update the xrange on the 2D plot to match the 1D"""
       self.plot_spec1D_current = False
       self.plot_spec2D_current = False
       self.plot_redshift_current = True
       
       self.mouse_x_redshift = self.plot_redshift.mapToView(pos).x()
       self.mouse_y_redshift = self.plot_redshift.mapToView(pos).y()
   
       self.setTitle_redshift()
      
   def mouseMoved_spec2D(self, pos):
       """Keep track of where the mouse and update the xrange on the 2D plot to match the 1D"""
       self.plot_spec1D_current = False
       self.plot_spec2D_current = True
       self.plot_redshift_current = False
       
       self.mouse_x_spec2D = self.plot_spec2D.mapToView(pos).x()
       self.mouse_y_spec2D = self.plot_spec2D.mapToView(pos).y()
       
          
   
   def setTitle_1D(self):
      
      self.plot_spec1D.setTitle('{} {}/{}   {:0.2f}, {:.2E}'.format(self.id, self.row, self.nRows, self.mouse_x_spec1D, self.mouse_y_spec1D))
   
   def setTitle_redshift(self):
      
      self.plot_redshift.setTitle('{} {}/{}   {:0.4f}, {:.2E}'.format(self.id, self.row, self.nRows, self.mouse_x_redshift, self.mouse_y_redshift))
   
      
   def draw(self):
      
      # Clear plots
      self.plot_redshift.clear()
      self.plot_spec1D.clear()
      
      
      if self.redshifted == 1:
         
         self.plot_redshift.plot(self.redshifts['z'], self.redshifts['chi2_pdf'], 
                                 pen=pg.mkPen('w', width=2), clear=True)
         self.plot_redshift.addItem(pg.InfiniteLine(self.z,
                                  pen=pg.mkPen('r', width=2, style=QtCore.Qt.DotLine)))
         self.plot_redshift.autoRange()
      self.plot_spec1D.plot(self.wave, self.flux1D*self.spec['mask'],
                        pen=pg.mkPen('w', width=2), clear=True)
      self.plot_spec1D.plot(self.wave, self.error1D*self.spec['mask'],
                        pen=pg.mkPen('c', width=2))
      # self.plot_spec1D.plot(self.wave, self.spec['model']*self.spec['mask'],
                        # pen=pg.mkPen('r', width=2))
      self.plot_spec1D.plot(self.wave, self.spec['model'],pen=pg.mkPen('r', width=2))
      # self.plot_spec1D.plot(self.wave, np.median(self.flux1D)*self.spec['mask'],pen=pg.mkPen('g', width=4)) #For debugging
      if self.param['Show raw']:
         self.plot_spec1D.plot(self.wave, self.flux1Draw*self.spec['mask'], 
            pen=pg.mkPen('w', width=1,style=QtCore.Qt.DotLine))
      #self.plot_spec1D.setYRange(np.percentile(self.flux1D, [5]),
      #                         np.percentile(self.flux1D, [99.9]),
      #                         padding=0)
                               
      self.setTitle_1D()
      self.setTitle_redshift()
      if self.param['Show lines']:
         
         features = self.features
         observedWaves = features['wave']*(1 + self.z)
         features = features[((observedWaves > np.min(self.wave)) & (observedWaves < np.max(self.wave))) | (features['list'] == 'sky')]
         
         for feature in features:
            
            if feature['list'] == 'sky':
               self.plot_spec1D.addItem(pg.InfiniteLine(feature['wave'],
                                        pen=pg.mkPen('g', width=2, style=QtCore.Qt.DotLine),
                                        label='{} {:0.1f}'.format(feature['name'], feature['wave']),
                                        labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
            elif feature['list'] == 'quasar':
               
               self.plot_spec1D.addItem(pg.InfiniteLine(feature['wave']*(1 + self.z),
                                        pen=pg.mkPen('b', width=2, style=QtCore.Qt.DotLine),
                                        label='{} {:0.1f}'.format(feature['name'], feature['wave']),
                                        labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
                                        
            elif feature['list'] == 'absorption':
               
               self.plot_spec1D.addItem(pg.InfiniteLine(feature['wave']*(1 + self.z),
                                        pen=pg.mkPen('r', width=2, style=QtCore.Qt.DotLine),
                                        label='{} {:0.1f}'.format(feature['name'], feature['wave']),
                                        labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
                                        
            elif feature['list'] == 'qsoals':
               
               self.plot_spec1D.addItem(pg.InfiniteLine(feature['wave']*(1 + self.z),
                                        pen=pg.mkPen('r', width=2, style=QtCore.Qt.DotLine),
                                        label='{} {:0.1f}'.format(feature['name'], feature['wave']),
                                        labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
                                        
            elif feature['list'] == 'emission':
               
               self.plot_spec1D.addItem(pg.InfiniteLine(feature['wave']*(1 + self.z),
                                        pen=pg.mkPen('y', width=2, style=QtCore.Qt.DotLine),
                                        label='{} {:0.1f}'.format(feature['name'], feature['wave']),
                                        labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))                        
         
      # self.plot_spec2D_view.addItem(self.plot_spec2D)

      # 2D spectrum

      self.plot_spec2D.setImage(self.flux2D, xvals=self.wave, 
                                levels=self.plot_spec2D_hist.getLevels(),#np.percentile(self.flux2D, [5, 99.5]), 
                                border=pg.mkPen('w', width=2))
      
      if self.version=='carpy':
         try: #need to remove those lines if they are there already...
            self.plot_spec2D_plot.removeItem(self.ext1)
            self.plot_spec2D_plot.removeItem(self.ext2)
            self.plot_spec2D_plot.removeItem(self.ext3)
         except:
            pass
         if self.param['Show trace']:
            try:
               pen=pg.mkPen('r',width=2,style=QtCore.Qt.DashLine)
               self.ext1=self.plot_spec2D_plot.plot(self.wave,self.extpos+self.exttrace,pen=pen)
               self.ext2=self.plot_spec2D_plot.plot(self.wave,self.extpos+self.exttrace+self.extaper/2,pen=pen)
               self.ext3=self.plot_spec2D_plot.plot(self.wave,self.extpos+self.exttrace-self.extaper/2,pen=pen)
            except:
               print('Trace file not loaded')

      self.objectsTable.selectRow(self.row-1)
      self.objectsTable.scrollToItem(self.objectsTable.item(self.row-1,0),QtGui.QAbstractItemView.PositionAtCenter)
      # temp=self.objectsTable.indexFromItem(self.objectsTable.item(self.row-1,0))
      # import pdb
      # pdb.set_trace()
      # self.objectsTable.scrollTo(temp,QtGui.QAbstractItemView.PositionAtCenter)
      



# Set up the command line argument parser
parser = argparse.ArgumentParser(description='Verify foreground/background quasars and measure absorbing gas around the foreground in the background quasar spectrum')
parser.add_argument('-m', metavar='mask name', type=str, help='name of the mask to be redshifted', required=True)
parser.add_argument('-xsize', metavar='xsize', type=int, help='xsize in pixels', default=2500)
parser.add_argument('-ysize', metavar='ysize', type=int, help='ysize in pixels', default=1500)
parser.add_argument('-v', metavar='version', type=str, help='input version (cosmos/carpy)', default='carpy')
args = parser.parse_args()

# Check for the 1d spectrum files.
if not os.path.isdir('{}_spec1D'.format(args.m)):
   print('Creating 1D spectrum files')
   createSpec1Dfiles(args.m,args.v)
     
else:
   print('1D spectrum files already present. Copying objects file to backup')
   shutil.copy('{}_spec1D/{}_objects.fits'.format(args.m,args.m),'{}_spec1D/{}_objects_bkp.fits'.format(args.m,args.m))
      
redshiftgui = ldss3_redshiftgui(args.m, args.xsize, args.ysize, args.v)