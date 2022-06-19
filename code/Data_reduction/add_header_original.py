#
# Copyright 2019-2022 Universidad Complutense de Madrid
#
# This file software has been employed to reduce infrared 
# raw data from the CIRCE camera at GTC (see Cabello et al. 2022)
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

from astropy.io import fits
import numpy as np


def read_raw_circe(infile1, infile2, outfile, debug=False):

    print("\n* Working with image:\n" + infile1)

    # read infile1
    hdulist = fits.open(infile1)
    main_header = hdulist[0].header
    # remove header keywords that are not FITS standard
    for wrong_keywords in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']:
        main_header.remove(wrong_keywords)
    # set keywords to be requested by imcombine
    main_header.set('TELESCOP', 'GTC', 'Telescope name (for imcombine)')
    airmass1 = main_header.get('AIRMASS1')
    airmass2 = main_header.get('AIRMASS2')
    airmass = (airmass1 + airmass2) / 2
    main_header.set('AIRMASS', airmass,
                    'Airmass (average of AIRMASS1 and AIRMASS2)',
                    before='AIRMASS1')
    hdulist.close()

    # read infile2
    hdulist = fits.open(infile2)
    result = hdulist[0].data
    hdulist.close()

    # save output file
    hdu = fits.PrimaryHDU(result, main_header)  # add main_header
    hdu.writeto(outfile, clobber=True)
    print("\nGenerating output file:\n" + outfile)


if __name__ == "__main__":

    list_object = """
    OB0001/object/0000885811-20160520-CIRCE-BROADBAND_IMAGE.fits
    OB0001/object/0000885812-20160520-CIRCE-BROADBAND_IMAGE.fits
    OB0001/object/0000885813-20160520-CIRCE-BROADBAND_IMAGE.fits
    OB0001/object/0000885814-20160520-CIRCE-BROADBAND_IMAGE.fits
    OB0001/object/0000885815-20160520-CIRCE-BROADBAND_IMAGE.fits
    OB0001/object/0000885816-20160520-CIRCE-BROADBAND_IMAGE.fits
    OB0001/object/0000885817-20160520-CIRCE-BROADBAND_IMAGE.fits
    OB0001/object/0000885818-20160520-CIRCE-BROADBAND_IMAGE.fits
    OB0001/object/0000885819-20160520-CIRCE-BROADBAND_IMAGE.fits
    OB0001/object/0000885820-20160520-CIRCE-BROADBAND_IMAGE.fits
    OB0001/object/0000885821-20160520-CIRCE-BROADBAND_IMAGE.fits
    """
    list_object = list_object.split()
    for infile1 in list_object:
        infile2 = infile1[:-5] + "_median_uc.fits"
        outfile = infile2[:-5] + "h.fits"
        read_raw_circe(infile1, infile2, outfile, debug=True)
