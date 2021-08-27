#!/usr/bin/env python
import sys
import numpy as np
from astropy.table import Table
import argparse
import redshift

# Set up the command line argument parser
parser = argparse.ArgumentParser(description='Compare HWC and SDJ MUSE redshifts')
parser.add_argument('-hwc', metavar='hwc', type=str, help='HWC catalog filename', required=True)
parser.add_argument('-sdj', metavar='sdj', type=str, help='SDJ catalog filename', required=True)

args = parser.parse_args()

hwc = Table.read(args.hwc, format='ascii')

print(hwc)

objects = Table.read('{}_spec1D/{}_objects.fits'.format(args.sdj, args.sdj))

objects.rename_column('class', 'class_SDJ')
objects.rename_column('quality', 'quality_SDJ')
objects.rename_column('redshift', 'redshift_SDJ')
objects.rename_column('comment', 'comment_SDJ')



objects['id_HWC'] = hwc['ID']


index = np.where(objects['id_HWC'] != objects['id'])[0]
if len(index) > 0:
   print('IDs do not match for all objects, check catalogs and try again')
   sys.exit()

print('All IDs match')

objects['rmag_auto'] = hwc['rmag_auto']
objects['drmag_auto'] = hwc['drmag_auto']

objects['redshift_HWC'] = hwc['redshift']
objects['quality_HWC'] = 0
index = np.where(objects['redshift_HWC'] > -1)
objects['quality_HWC'][index] = 2

index = np.where(objects['quality_SDJ'] <= 0)[0]
objects['redshift_SDJ'][index] = -2.0




objects['dv'] = redshift.dv(objects['redshift_SDJ'], objects['redshift_HWC'])


objects['redshift_SDJ'].info.format = '0.4f'
objects['dv'].info.format = '0.0f'



index = np.where((objects['quality_SDJ'] == 0) & (objects['quality_HWC'] == 0))
objects['dv'][index] = 0.0


objects['check'] = 0

index = np.where((objects['quality_SDJ'] != objects['quality_HWC']) | (np.abs(objects['dv']) > 20.0))
objects['check'][index] = 1

index = np.where((objects['quality_SDJ'] == 2) & (objects['class_SDJ'] == 'star') & (objects['redshift_HWC'] == 0.0))
objects['check'][index] = 0


print(objects)


objects = objects['row', 'id', 'name', 'ra', 'dec', 'radius', 'rmag_auto', 'drmag_auto',
                  'check', 'quality_SDJ', 'quality_HWC', 'class_SDJ', 'redshift_SDJ', 'redshift_HWC', 'dv', 'comment_SDJ']
print(objects)


objects_check = objects[objects['check'] == 1]
print(objects_check)

objects.write(args.hwc.replace('.zcat', '_comparison.dat'), overwrite=True, format='ascii.fixed_width', delimiter='')
