#!/usr/bin/env python
import sys
import numpy as np
from astropy.table import Table, hstack, vstack, join
from matplotlib import pyplot as p
import argparse

c_kms = 299792.458

def readMask(filename, number):
   
   
   objects = Table.read(filename)
   #objects.remove_columns(['extpos', 'extflag', 'alignbox', 'extaper'])
   objects['redshift'].info.format = '0.5f'
   
   objects.rename_column('class', 'class{}'.format(number))
   objects.rename_column('redshift', 'redshift{}'.format(number))
   objects.rename_column('quality', 'quality{}'.format(number))
   objects.rename_column('comment', 'comment{}'.format(number))
   
   return objects


# Set up the command line argument parser
parser = argparse.ArgumentParser(description='Compare two redshift catalogs assigned by different people for the same mask and flag objects with differences')
parser.add_argument('-o1', metavar='obect catalog 1', type=str, help='mask object file from user 1', required=True)
parser.add_argument('-o2', metavar='obect catalog 2', type=str, help='mask object file from user 2', required=True)

args = parser.parse_args()

objects1 = readMask(args.o1, 'A')
objects2 = readMask(args.o2, 'B')

print(objects1)
print(objects2)

index = np.where(objects1['id'] == objects2['id'])[0]


if (len(index) == len(objects1)) & (len(index) == len(objects2)):
   
   print('All lengths and object IDs match, proceeding')
   
   
   print(objects1)
   print(objects2)
   
   #objects2.remove_columns(['row', 'id'])
   objects = join(objects1, objects2, keys=['row', 'id'])
   objects['dv(A-B)'] = (objects['redshiftA'] - objects['redshiftB'])/(1 + objects['redshiftA'])*c_kms
   objects = objects['row', 'id', 'qualityA', 'qualityB', 'redshiftA', 'redshiftB', 'dv(A-B)', 'classA', 'classB', 'commentA', 'commentB']
   print(objects)
   
   objects = objects[(objects['qualityA'] >= 1) | (objects['qualityB'] >= 1)]
   
   objects_bigdisagreement = objects[(objects['qualityA'] != objects['qualityB']) | (objects['classA'] != objects['classB']) | (np.abs(objects['dv(A-B)']) > 200.0)]
   
   objects_minordisagreement = objects[(objects['qualityA'] == objects['qualityB']) & (objects['classA'] == objects['classB']) & (np.abs(objects['dv(A-B)']) > 60.0) & (np.abs(objects['dv(A-B)']) <= 200.0)]
   
   objects_disagreement = vstack([objects_bigdisagreement, objects_minordisagreement])
   
   
   
   print('Objects with significant disagreements')
   print(objects_bigdisagreement)
   
   print('Objects with minor disagreements')
   print(objects_minordisagreement)
   
   print('Objects with disagreements')
   print(objects_disagreement)
   
   filename1 = args.o1
   savename = filename1.split('/')[-1]
   savename = savename.replace('_objects.fits', '_z_comparisons.txt')
   print(savename)
   
   objects_disagreement.write(savename, overwrite=True, format='ascii.fixed_width', delimiter='')
   
else:
   
   print('Table lengths and/or IDs do not match')
   sys.exit()



