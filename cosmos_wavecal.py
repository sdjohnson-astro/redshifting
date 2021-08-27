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
import lmfit
from scipy import stats


def getpolyfit(wave, poly):
   
   deg = len(poly) - 1
   i = 0
   polyfit = np.zeros(len(wave))
   for i in range(len(poly)):
      
      print(i, deg - i)
      
      polyfit = polyfit + poly[i]*wave**(deg - i)
      
   return polyfit


def arcfit(wave, a, b, wave_fit, amplitude, sigma):
         
      gaussian = amplitude*np.exp(-(wave - wave_fit)**2/2/sigma**2)
      continuum = a*wave + b
         
      return gaussian + continuum

def getspec1Dname(mask, row, id):
   
   return '{}_spec1D/{}_{}_{}.fits'.format(mask, mask, row, id)


class cosmos_wavecalgui:
   """Defining the GUI class for cosmos wavelength calibration updates"""
   
   def __init__(self, mask, xsize=1000, ysize=1000, version='carpy'):
      
      self.mask = mask
      self.xsize = xsize
      self.ysize = ysize
      self.dW = 15.0
      self.degree = 4
      
      self.lines = Table.read('HeNeAr_check.dat', format='ascii')
      
      lines_fit = np.copy(self.lines)
      lines_fit = Table(lines_fit)
      lines_fit['a'] = 0.0
      lines_fit['b'] = 0.0
      lines_fit['wave_fit'] = 0.0
      lines_fit['amplitude'] = 0.0
      lines_fit['sigma'] = 0.0
      lines_fit['dW'] = 0.0
      self.lines_fit = lines_fit
      
      self.path = '{}_spec1D'.format(mask)
      self.objects = Table.read('{}/{}_objects.fits'.format(self.path, self.mask))
      print(self.objects)
      
      # Set the initial row number to zero
      self.row = 1
      self.nRows = len(self.objects)
      self.smoothing = 1
      
      
      # Get the GUI ready
      self.app = QtGui.QApplication([])       # Always start by initializing Qt
      self.widget = QtGui.QWidget()       # Define a top-level widget
      
      # Set the widget size
      self.widget.resize(self.xsize, self.ysize)
      
      
      # Set the background plotting widget
      
      self.plot_solution = pg.PlotWidget()
      self.plot_solution.getAxis('bottom').setPen(pg.mkPen('w', width=2))
      self.plot_solution.getAxis('top').setPen(pg.mkPen('w', width=2))
      self.plot_solution.getAxis('left').setPen(pg.mkPen('w', width=2))
      self.plot_solution.getAxis('right').setPen(pg.mkPen('w', width=2))
      self.plot_solution.getAxis('bottom').setStyle(tickLength=-15)
      self.plot_solution.getAxis('top').setStyle(tickLength=-15)
      self.plot_solution.getAxis('left').setStyle(tickLength=-15)
      self.plot_solution.getAxis('right').setStyle(tickLength=-15)
      self.plot_solution.showAxis('right')
      self.plot_solution.showAxis('top')
      self.plot_solution.setLabel('bottom', 'wavelength')
      self.plot_solution.setLabel('left', 'dW [&#8491;]')


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
      self.plot_spec1D.setLabel('bottom', 'arc flux')

      self.plot_spec1D.setXLink(self.plot_solution)

      self.mouse_x_spec1D = 0.0
      self.mouse_y_spec1D = 0.0
      self.mouse_x_solution = 0.0
      self.mouse_y_solution = 0.0
      
      # Setup the layout
      self.layout = QtGui.QGridLayout()
      self.widget.setLayout(self.layout)
      
      # Add listeners
      self.plot_spec1D.scene().sigMouseMoved.connect(self.mouseMoved_spec1D)
      self.plot_spec1D.keyPressEvent = self.keypress_spec1D
      
      self.plot_solution.scene().sigMouseMoved.connect(self.mouseMoved_solution)
      self.plot_solution.keyPressEvent = self.keypress_solution
      
      self.objectsTable = pg.TableWidget(editable=False, sortable=False)
      self.objectsTable.setFormat('%0.5f', 2)
      self.setTable()
      
      self.objectsTable.setSelectionBehavior(QtGui.QAbstractItemView.SelectRows)
      self.objectsTable.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
      self.objectsTable.doubleClicked.connect(self.goToObject)
      
      self.linesTable_fit = pg.TableWidget(editable=False, sortable=False)
      self.linesTable_fit.setFormat('%0.2f', 0)
      self.linesTable_fit.setFormat('%0.2f', 1)
      self.linesTable_fit.setFormat('%0.2f', 2)
      self.linesTable_fit.setFormat('%0.2f', 3)
      self.linesTable_fit.setFormat('%0.2f', 4)
      self.linesTable_fit.setFormat('%0.2f', 5)
      self.linesTable_fit.setFormat('%0.2f', 6)
      
      self.setTable_lines()
      
      
      
      # Add plot_bg to the layout
      self.layout.addWidget(self.plot_solution, 0, 0)
      self.layout.addWidget(self.linesTable_fit, 0, 1)
      self.layout.addWidget(self.objectsTable, 1, 1)
      self.layout.addWidget(self.plot_spec1D, 1, 0)
      #self.layout.addWidget(self.comment_text, 4, 0)
      
      
      self.layout.setColumnStretch(0, 0)
      self.layout.setColumnStretch(0, 1)
      self.layout.setColumnStretch(1, 0)
      #self.layout.setColumnStretch(1, 1)
      
      self.layout.setRowStretch(0, 0)
      self.layout.setRowStretch(0, 1)
      self.layout.setRowStretch(1, 0)
      self.layout.setRowStretch(1, 1)
      
      self.setSpec()
      
      self.fitAll()
      self.fitResiduals()
      dW = self.lines_fit['dW'] - getpolyfit(self.lines_fit['wave'], self.poly)
      sigma_robust = stats.median_absolute_deviation(dW)
      print('robust_sigma = {}'.format(sigma_robust))
      index = np.where(np.abs(dW) > 5.0*sigma_robust)
      self.lines_fit['a'][index]         = 0.0
      self.lines_fit['b'][index]         = 0.0
      self.lines_fit['wave_fit'][index]  = 0.0
      self.lines_fit['amplitude'][index] = 0.0
      self.lines_fit['sigma'][index]     = 0.0
      self.fitResiduals()
      
      self.draw()
      
      self.widget.show()
      self.app.exec_()
      
      
      
   def setTable(self):
   
      self.objectsTable.setData(np.array(self.objects['row', 'id', 'class', 'redshift', 'quality']))

   def setTable_lines(self):
   
      self.linesTable_fit.setData(np.array(self.lines_fit))

   
   def goToObject(self):
   
      #print('Going to object...')
      self.row = self.objectsTable.selectedItems()[0].row()+1
      self.setSpec()
      # self.draw()
   
   
   def keypress_solution(self, event):
   
      if self.plot_solution_current: 
      
         #print('')
         #print(event.text())
         #print(event.key())
         #print('{:0.4f}, {:0.2f}'.format(self.mouse_x_solution, self.mouse_y_solution))
   
         if event.text() in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
            
            self.degree = int(event.text())
            self.setTitle_solution()
            
         if event.text() in ['t','T']:

            self.fitResiduals()
            
         if event.text() in ['r','R']:

            self.redoWavecal()
         
   
         if (event.text() == 'd'):
            
            self.deleteNearest()
   
         if (event.text() == 'n'):
         
            self.advance(1)
      
         if (event.text() == 'N'):
         
            self.advance(10)
         
         if (event.text() == 'b'):
         
            self.advance(-1)
         
         if (event.text() == 'B'):
         
            self.advance(-10)
      
         if event.text() == '[':
            self.panx_solution(-1.0/3.0)
         
         if event.text() == ']':
            self.panx_solution(1.0/3.0)
         
         if event.text() == '{':
            self.panx_solution(-1.0)
         
         if event.text() == '}':
            self.panx_solution(1.0)
         
         if event.text() == 'x':
            self.zoomxy_solution(1.0/1.5, 1.0)
      
         if event.text() == 'X':
            self.zoomxy_solution(1.5, 1.0)
         
         if event.text() == 'y':
            self.zoomxy_solution(1.0, 1.0/1.5)
      
         if event.text() == 'Y':
            self.zoomxy_solution(1.0, 1.5)
         
         if (event.text() == 'w') | (event.text() == 'W'):
            self.plot_spec1D.autoRange()
   
   
      
         
   def zoomxy_solution(self, scalex, scaley):
      """Zoom in or out in wavelength (x) and/or flux (y)"""
   
      xRange = self.plot_solution.getViewBox().state['viewRange'][0]
      x0 = xRange[0]
      x1 = xRange[1]
      xRange = (x1 - x0)*scalex
      x0_new = self.mouse_x_solution - xRange/2.0
      x1_new = self.mouse_x_solution + xRange/2.0
   
      self.plot_solution.setXRange(x0_new, x1_new, padding=0)
   
      yRange = self.plot_solution.getViewBox().state['viewRange'][1]
      y0 = yRange[0]
      y1 = yRange[1]
      yRange = (y1 - y0)*scaley
      y0_new = self.mouse_y_solution - yRange/2.0
      y1_new = self.mouse_y_solution + yRange/2.0
   
      self.plot_solution.setYRange(y0_new, y1_new, padding=0)

   def zoom_default(self):
      q1,q2=np.percentile(self.arc1D,[0.1,99.9])
      q1 = np.min([0.0, q1])
      self.plot_spec1D.setYRange(-0.05*q2,1.1*q2,padding=0)
      self.plot_spec1D.setXRange(np.min(self.wave), np.max(self.wave), padding=0)

   
   # def zoom_solution(self, key, ):
   #    """Zoom in or out in wavelength (x) and/or flux (y)"""
   
   #    xRange = self.plot_solution.getViewBox().state['viewRange'][0]
   #    x0 = xRange[0]
   #    x1 = xRange[1]
   #    xRange = (x1 - x0)*scalex
   #    x0_new = self.mouse_x_solution - xRange/2.0
   #    x1_new = self.mouse_x_solution + xRange/2.0
   
   #    self.plot_solution.setXRange(x0_new, x1_new, padding=0)
   
   #    yRange = self.plot_solution.getViewBox().state['viewRange'][1]
   #    y0 = yRange[0]
   #    y1 = yRange[1]
   #    yRange = (y1 - y0)*scaley
   #    y0_new = self.mouse_y_solution - yRange/2.0
   #    y1_new = self.mouse_y_solution + yRange/2.0
   
   #    self.plot_solution.setYRange(y0_new, y1_new, padding=0)
   
   def panx_solution(self, scalex):
      """Pan in the wavelength direction"""
   
      xRange = self.plot_solution.getViewBox().state['viewRange'][0]
      x0 = xRange[0]
      x1 = xRange[1]
      shift = scalex*(x1 - x0)
      x0_new = x0 + shift
      x1_new = x1 + shift

      self.plot_solution.setXRange(x0_new, x1_new, padding=0)
      
      
      
      
   def keypress_spec1D(self, event):
   
   
      if self.plot_spec1D_current:
   
         #print(event.text())
         #print(event.key())
         #print('')
         
         if event.text() in ['r','R']:

            self.redoWavecal()
         
         if event.text() in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
            
            self.degree = int(event.text())
            self.setTitle_solution()
   
         if (event.text() == '+') | (event.text() == '='):
            
            self.dW = self.dW + 1
            self.setTitle_1D()
            
         if (event.text() == '-') | (event.text() == '_'):
            
            self.dW = self.dW - 1
            self.setTitle_1D()
   
         if (event.text() == 'f'):
            
            self.fitNearest()
            
         if (event.text() == 'F'):
            
            self.fitAll()
            
         if (event.text() == 'd'):
            
            self.deleteNearest()
   
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

         if event.text() in ['t','T']:

            self.fitResiduals()
      
         if (event.text() == 'w') | (event.text() == 'W'):
            self.zoom_default()
            # self.plot_spec1D.autoRange()
            # self.updateXrange_1D()            
      
      
         if event.text() == 'm':
      
            self.changeMask(0)
      
         if event.text() == 'M':
      
            # self.changeMask(0)
            self.autoMask()
      
         if event.text() == 'u':
      
            self.changeMask(1)
      
         if event.text() == 'U':
      
            self.changeMask(1)
      
      
      
         # if return or enter is pressed then save.
         if (event.key() == 16777220) | (event.key() == 16777221) | (event.key() == 96) | (event.key() == 126):
      
            self.save()


   def deleteNearest(self):
      

      if self.plot_spec1D_current:
         print('Deleting from spec1D')
         dW = np.abs(self.lines_fit['wave'] - self.mouse_x_spec1D)
      elif self.plot_solution_current:
         print('Deleting from solution')
         dW = np.abs(self.lines_fit['wave'] - self.mouse_x_solution)
      
      minIndex = np.argmin(dW)
      wave_closest = self.lines[minIndex]['wave']
      #print('wave_closest = {}'.format(wave_closest))
      
      self.lines_fit['a'][minIndex]         = 0.0
      self.lines_fit['b'][minIndex]         = 0.0
      self.lines_fit['wave_fit'][minIndex]  = 0.0
      self.lines_fit['amplitude'][minIndex] = 0.0
      self.lines_fit['sigma'][minIndex]     = 0.0
      
      
      self.setTable_lines()
      
      self.draw()

   def redoWavecal(self):
      
      self.wave = self.wave - self.polyfit
      self.spec['wave'] = self.wave
      self.polyfit = np.zeros(len(self.wave))
      self.poly = [0.0]
      self.draw()

   def fitResiduals(self):
      
      thisLines = self.lines_fit
      thisLines = thisLines[thisLines['wave_fit'] > 0.0]
      self.poly = np.polyfit(thisLines['wave'], thisLines['dW'], deg=self.degree)
      print('Polynomial fit performed')
      print(thisLines)
      
      print(self.poly)
      self.polyfit = getpolyfit(self.wave, self.poly)
      self.draw()

   def fitNearest(self):
      
      dW = np.abs(self.lines['wave'] - self.mouse_x_spec1D)
      
      minIndex = np.argmin(dW)
      wave_closest = self.lines[minIndex]['wave']
      #print('wave_closest = {}'.format(wave_closest))
      
      self.fitArcLine(wave_closest)
      
      
   def fitAll(self):
      
      print('Fitting all lines')
      for line in self.lines:
         
         
         
         try:
            self.fitArcLine(line['wave'], draw=False)
         except:
            print('Fit unsuccessful for {:0.2f}'.format(line['wave']))
      
      index = np.where(np.abs(self.lines_fit['dW']) > 10.0)
      self.lines_fit['a'][index] = 0.0
      self.lines_fit['b'][index] = 0.0
      self.lines_fit['wave_fit'][index] = 0.0
      self.lines_fit['amplitude'][index] = 0.0
      self.lines_fit['sigma'][index] = 0.0
      self.lines_fit['dW'][index]  = 0.0
      
      self.fitResiduals()
      
      self.draw()
      
      
      
   def fitArcLine(self, wave_closest, draw=True):
      
      
      dW = self.dW
      thisWave = np.array(self.wave)
      thisArc1D = np.array(self.arc1D)
      
      index = np.where(np.abs(thisWave - wave_closest) <= dW)[0]
      thisWave = thisWave[index]
      thisArc1D = thisArc1D[index]
      
      if np.max(thisArc1D) > 0.0:
      
         parameters = lmfit.Parameters()
         parameters.add_many(('a',         np.min(thisArc1D),                 True,  None,  None,  None),
                             ('b',         0.0,                 True,  None,  None,  None),
                             ('wave_fit',  wave_closest,        True,  wave_closest-dW, wave_closest+dW,  None),
                             ('amplitude', np.max(thisArc1D) - np.min(thisArc1D),   False, None,  None,   None),
                             ('sigma',     5.0,                 True,  1.00,  10.0,   None))
         
         arc_model = lmfit.Model(arcfit, missing='drop')
         result = arc_model.fit(thisArc1D,
                                 wave=thisWave,
                                 params=parameters)
         
         #print('success = {}'.format(result.success))
         #print(result.success)
         #print(result.fit_report())
         if result.success:
            index = np.where(np.abs(self.lines_fit['wave'] - wave_closest) <= dW)[0]
            self.lines_fit['a'][index] = result.best_values['a']
            self.lines_fit['b'][index] = result.best_values['b']
            self.lines_fit['wave_fit'][index] = result.best_values['wave_fit']
            self.lines_fit['amplitude'][index] = result.best_values['amplitude']
            self.lines_fit['sigma'][index] = result.best_values['sigma']
            self.lines_fit['dW'][index]  = self.lines_fit['wave_fit'][index]  - self.lines_fit['wave'][index] 
            
            self.setTable_lines()
         else:
            print('fit not successful for {:0.2f}'.format(wave_closest))
      else:
         print('No valid data for this arc line')
      
      if draw:
         self.draw()
      
      
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
      #if self.redshifted == 1:
      #   savename = getredshift1Dname(self.mask, self.objects[self.row-1]['row'],
      #                            self.objects[self.row-1]['id'])
      #   self.redshifts.write(savename, overwrite=True)
   
      print('Saved')   

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
   
      


#   def trimxy_spec(self, key):
#      """Trim plotting region in wavelength (x) and/or flux (y)"""
#      if key in 'eE':
#         xRange = self.plot_spec1D.getViewBox().state['viewRange'][0]
#         x0 = xRange[0]
#         x1 = xRange[1]
#         xnew = self.mouse_x_spec1D
#         if key=='E': x1=xnew
#         else: x0=xnew
#         self.plot_spec1D.setXRange(x0, x1, padding=0)
#      elif key in 'tT':
#         yRange = self.plot_spec1D.getViewBox().state['viewRange'][1]
#         y0 = yRange[0]
#         y1 = yRange[1]
#         ynew = self.mouse_y_spec1D
#         if key=='t': y1=ynew
#         else: y0=ynew
#         self.plot_spec1D.setYRange(y0, y1, padding=0)
#      #self.updateXrange_1D()

      
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
   
      #self.updateXrange_1D()

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
      #self.updateXrange_1D()



   def mouseMoved_spec1D(self, pos):
       """Keep track of where the mouse and update the xrange on the 2D plot to match the 1D"""
       self.plot_spec1D_current = True
       self.plot_solution_current = False
    
       self.mouse_x_spec1D = self.plot_spec1D.mapToView(pos).x()
       self.mouse_y_spec1D = self.plot_spec1D.mapToView(pos).y()
    
       #self.updateXrange_1D()
    
       self.setTitle_1D()
    
 

   def advance(self, delt):
   
      self.save()
      self.row = self.row + delt
      print(self.row)
      if self.row < 1:
         self.row = self.nRows
      if self.row > self.nRows:
         self.row = 1
      self.setSpec()
      
      self.fitAll()
      self.fitResiduals()
      
      dW = self.lines_fit['dW'] - getpolyfit(self.lines_fit['wave'], self.poly)
      sigma_robust = stats.median_absolute_deviation(dW)
      print('robust_sigma = {}'.format(sigma_robust))
      index = np.where(np.abs(dW) > 5.0*sigma_robust)
      self.lines_fit['a'][index]         = 0.0
      self.lines_fit['b'][index]         = 0.0
      self.lines_fit['wave_fit'][index]  = 0.0
      self.lines_fit['amplitude'][index] = 0.0
      self.lines_fit['sigma'][index]     = 0.0
      self.fitResiduals()
         
      #print('{}/{}'.format(self.row, self.nRows))
      self.draw()
   
   def setSpec(self, autoRange=True):
      
      print('Need to setup spectrum')
      """Set the spectrum to current row"""
      #self.z = self.objects[self.row-1]['redshift']
      #self.comment_text.setText(self.objects[self.row-1]['comment'])
      #
      self.id=self.objects[self.row-1]['id']
      #   if not os.path.isfile(getspec1Dname(self.mask, self.row, self.id)):
      #      print('Files for {} do not exist. Moving on.'.format(self.id))
      #      self.redshifted=0
      #      self.advance(1)
      #      return
      #   self.flux2D=fits.getdata(getspec2Dname(self.mask,self.id)).transpose()
      #
      #elif self.version=='cosmos':
      #   # Get the apnum header parameter
      #   self.apnum1D = self.header1D['APNUM{}'.format(self.row)]
      #   self.apnum2D = self.header2D['APNUM{}'.format(self.row)]
      #   self.apnum1Darray = self.apnum1D.split(' ')
      #   self.apnum2Darray = self.apnum2D.split(' ')
      #   self.id = self.apnum1Darray[1]
      #   self.y0=int(self.header2D['CSECT{}A'.format(self.row)])
      #   self.y1=int(self.header2D['CSECT{}B'.format(self.row)])
      #   self.y0 = int(float(self.apnum2Darray[2]))
      #   self.y1 = int(float(self.apnum2Darray[3]))
      #   self.flux2D = self.spec2Darray[self.y0:self.y1, :].transpose()
      #
      self.spec = fits.getdata(getspec1Dname(self.mask, self.objects[self.row-1]['row'], self.id))
      #
      self.wave = self.spec['wave']
      self.flux1D = self.spec['flux']
      self.flux1Draw = self.spec['raw']
      self.error1D = self.spec['error']
      self.error1Draw=self.spec['rawerr']
      self.model1D = self.spec['model']
      self.flat1D = self.spec['flat']
      self.arc1D = self.spec['arc']
      
      self.poly = [0.0]
      self.polyfit = np.zeros(len(self.wave))
      
      self.lines_fit['a'] = 0.0
      self.lines_fit['b'] = 0.0
      self.lines_fit['wave_fit'] = 0.0
      self.lines_fit['amplitude'] = 0.0
      self.lines_fit['sigma'] = 0.0
      self.lines_fit['dW'] = 0.0
      
      self.poly = [0.0]
      self.polyfit = np.zeros(len(self.wave))
      
      #
      #self.extpos = self.objects[self.row-1]['extpos']
      #self.extaper = self.objects[self.row-1]['extaper']
      #self.smoothSpec()
      #
      ## Check for redshift filename and read in if present.
      #redshiftFilename = getredshift1Dname(self.mask, self.objects[self.row-1]['row'],
      #                                     self.objects[self.row-1]['id'])
      #if os.path.isfile(redshiftFilename):
      #   self.redshifted = 1
      #   self.redshifts = Table.read(redshiftFilename)
      #   #print('Already redshifted')
      #
      #else:
      #   self.redshifted = 0
      #   self.redshifts = None
      #
      self.draw()
      #self.plot_solution.autoRange()
      if autoRange:
         self.plot_spec1D.autoRange()
         self.zoom_default()
   

   def draw(self):
      
      # Clear plots
      self.plot_solution.clear()
      self.plot_spec1D.clear()
      
      
      self.plot_spec1D.plot(self.wave, self.arc1D,
                        pen=pg.mkPen('w', width=2), clear=True)
       
      
                                
      self.setTitle_1D()
      self.setTitle_solution()
 
      self.objectsTable.selectRow(self.row-1)
      self.objectsTable.scrollToItem(self.objectsTable.item(self.row-1,0),QtGui.QAbstractItemView.PositionAtCenter)
      # temp=self.objectsTable.indexFromItem(self.objectsTable.item(self.row-1,0))
      # import pdb
      # pdb.set_trace()
      # self.objectsTable.scrollTo(temp,QtGui.QAbstractItemView.PositionAtCenter)

      for line in self.lines:
         
         self.plot_spec1D.addItem(pg.InfiniteLine(line['wave'],
                                  pen=pg.mkPen('y', width=2, style=QtCore.Qt.DotLine),
                                  label='{:0.1f}'.format(line['wave']),
                                  labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
      
      
                                  
      for line in self.lines_fit:
         
         dW = self.dW
         if line['wave_fit'] > 0.0:
            
            self.plot_spec1D.addItem(pg.InfiniteLine(line['wave'],
                                     pen=pg.mkPen('g', width=2, style=QtCore.Qt.DotLine),
                                     label='{:0.1f}'.format(line['wave']),
                                     labelOpts={'position':0.8, 'rotateAxis':[1, 0]}))
            
            thisWave = self.wave
            thisWave = thisWave[(thisWave >= line['wave'] - dW) & (thisWave <= line['wave'] + dW)]
            thisModel = arcfit(thisWave, line['a'], line['b'], line['wave_fit'],
                               line['amplitude'], line['sigma'])
            self.plot_spec1D.plot(thisWave, thisModel,
                              pen=pg.mkPen('r', width=2))
                              
                              
      self.plot_solution.addItem(pg.InfiniteLine(0.0, angle=0.0,
                               pen=pg.mkPen('w', width=1, style=QtCore.Qt.DotLine)))
      
      self.plot_solution.plot(self.wave, self.polyfit,
                        pen=pg.mkPen('r', width=2), clear=True)
                              
      thisLines = self.lines_fit
      index = np.where(thisLines['wave_fit'] > 0)[0]
      if len(index) > 0:
         print('Plotting wavecal residuals')
         thisLines = thisLines[index]
         
         self.plot_solution.plot(thisLines['wave'], thisLines['dW'],
                           pen=None, symbol='o')
         
         #self.plot_solution.autoRange()
   
   def mouseMoved_solution(self, pos):
       """Keep track of where the mouse and update the xrange on the 2D plot to match the 1D"""
       self.plot_spec1D_current = False
       self.plot_spec2D_current = False
       self.plot_solution_current = True
    
       self.mouse_x_solution = self.plot_solution.mapToView(pos).x()
       self.mouse_y_solution = self.plot_solution.mapToView(pos).y()

       self.setTitle_solution()
   
     
    

   def setTitle_1D(self):
   
      self.plot_spec1D.setTitle('{} {}/{}   {:0.2f}, {:.2E}   dW_fit={:0.2f}'.format(self.id, self.row, self.nRows, self.mouse_x_spec1D, self.mouse_y_spec1D, self.dW))

   def setTitle_solution(self):
   
      self.plot_solution.setTitle('{} {}/{}   {:0.4f}, {:.2E}   order={}'.format(self.id, self.row, self.nRows, self.mouse_x_solution, self.mouse_y_solution, self.degree))

      
      

# Set up the command line argument parser
parser = argparse.ArgumentParser(description='Verify and improve COSMOS wavelength calibrations')
parser.add_argument('-m', metavar='mask name', type=str, help='name of the mask to be checked', required=True)
parser.add_argument('-xsize', metavar='xsize', type=int, help='xsize in pixels', default=2500)
parser.add_argument('-ysize', metavar='ysize', type=int, help='ysize in pixels', default=1500)
args = parser.parse_args()

# Check for the 1d spectrum files.
if not os.path.isdir('{}_spec1D'.format(args.m)):
   print('1D spectrum files not yet present. Run cubs_solutioning GUI before trying to fix wavelength calibrations')
   sys.exit()
     
else:
   print('1D spectrum files already present.')
      
wavecalgui = cosmos_wavecalgui(args.m, args.xsize, args.ysize)