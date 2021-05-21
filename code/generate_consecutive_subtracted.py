#
# Copyright 2019-2021 Universidad Complutense de Madrid
#
# This file software has been employed to reduce infrared 
# raw data from the CIRCE camera at GTC
#
# Authors: Cristina Cabello (criscabe@ucm.es), 
#          Nicolás Cardiel (cardiel@ucm.es)
#          Jesús Gallego (j.gallego@ucm.es)
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from __future__ import division
from __future__ import print_function

import argparse
from astropy.io import fits
import numpy as np
import re


def readfits(infile, iext=0):
    """Read FITS file.

    Parameters
    ----------
    infile : string
        Input file name to be read.
    iext : int
        Extension number (0=first HDU)

    Returns
    -------
    image2d : 2d numpy array, floats
        Data array.
    main_header: fits header
        Primary header.
    header : fits header
        Header corresponding to the extension read. When the extension
        is zero, 'header' and 'main_header' are the same.

    """

    hdulist = fits.open(infile)
    main_header = hdulist[0].header
    header = hdulist[iext].header
    image2d = hdulist[iext].data.astype(np.float)
    hdulist.close()

    return image2d, main_header, header


if __name__ == "__main__":

    # define expected command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_list", 
                        help='"txt file with list of initial file names"')
    parser.add_argument("nramps",
                        help="number of ramps",type=int)
    args = parser.parse_args()

    nramps = int(args.nramps)
    
    # read list with file names
    list_files = np.genfromtxt(args.input_list, dtype=[('filename', '|S100')])
    filename = list_files['filename']

    nfiles = filename.size

    filename_cs = np.copy(filename)

    # subtract consecutive images
    for k in range(nfiles):
        if k < nfiles - nramps:
            print("image1: ", filename[k])
            print("image2: ", filename[k+nramps])
            image1, main_header, dum = readfits(filename[k])
            image2, dum, dum = readfits(filename[k + nramps])
        else:
            print("image1: ", filename[k])
            print("image2: ", filename[k-nramps])
            image1, main_header, dum = readfits(filename[k])
            image2, dum, dum = readfits(filename[k - nramps])
        image = image1 - image2
        
        outfile = re.sub('.fits$', '_cs.fits', filename[k])
        print("Generating:" + outfile)
        hdu = fits.PrimaryHDU(image, main_header)
        hdu.writeto(outfile, clobber=True)
        filename_cs[k] = outfile
   
    output_txt = re.sub('.txt$', '_cs.txt', args.input_list)
    np.savetxt(output_txt, filename_cs, fmt="%s") 
    print("\nExporting list with new files to:", output_txt)
