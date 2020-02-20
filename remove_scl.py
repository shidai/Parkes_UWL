import math
import numpy as np

import matplotlib.pyplot as plt
import argparse

from astropy.io import fits

parser = argparse.ArgumentParser(description='Remove scales in psrfits file')
parser.add_argument('-f',  '--input_file',  metavar='Input file name',  nargs='+', required=True, help='Input file name')
parser.add_argument('-o',  '--output_file', metavar='Output file name', nargs='+', required=True, help='Output file name')

args = parser.parse_args()

in_file = args.input_file[0]
hdulist = fits.open(in_file)

############## update weights ###############
scl = hdulist['SUBINT'].data['DAT_SCL']
shape = scl.shape
scl_new = np.ones(shape)
hdulist['SUBINT'].data['DAT_SCL'] = scl_new

off = hdulist['SUBINT'].data['DAT_OFFS']
shape = off.shape
off_new = np.ones(shape)
hdulist['SUBINT'].data['DAT_OFFS'] = off_new
############## Output to a new file ###############
out_file = args.output_file[0]

hdulist.writeto(out_file, overwrite=True, checksum=True)

hdulist.close()
