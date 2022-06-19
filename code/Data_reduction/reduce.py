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

import argparse
from astropy.io import fits
import numpy as np
import re



def read_raw_circe(infile, nramps, outfile, debug=0):
    """Read raw CIRCE ramp and generate coadded result.

    Parameters
    ----------
    infile : string
        Input file name containing the individual reads along the
        different ramps.
    nramps : int
        Number of ramps.
    outfile : string
        Output file name corresponding to the result of coadding all
        the ramps.
    debug : int
        debug=0: no additional information is generated
        debug=1: save ramps
        debug=2: save ramps and individual reads in the ramps

    """

    print("\n* Working with image:\n" + infile)
    hdulist = fits.open(infile)
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
    lista_rampas = []
    for i in range(nramps):
        print("   - reading ramp " + str(i+1))
        image2d_1 = hdulist[2*i+1].data.astype(np.float)
        image2d_2 = hdulist[2*i+2].data.astype(np.float)
        if debug==2:
            hdu = fits.PrimaryHDU(image2d_1)
            dumfile = infile[:-5] + "_ext" + str(2*i+1) + ".fits"
            hdu.writeto(dumfile, clobber=True)
            hdu = fits.PrimaryHDU(image2d_2)
            dumfile = infile[:-5] + "_ext" + str(2*i+2) + ".fits"
            hdu.writeto(dumfile, clobber=True)
        result = image2d_2 - image2d_1
        if debug==1 or debug==2:
            hdu = fits.PrimaryHDU(result)
            dumfile = infile[:-5] + "_ramp" + str(i+1) + ".fits"
            hdu.writeto(dumfile, clobber=True)
        lista_rampas.append(result)
    hdulist.close()

    result = np.copy(lista_rampas[0])
    if nramps > 1:
        for i in range(1, nramps):
            result += lista_rampas[i]

    result /= nramps
    # save output file
    hdu = fits.PrimaryHDU(result, main_header)  # add main_header
    hdu.writeto(outfile, clobber=True)
    print("\nGenerating output file:\n" + outfile)

def average_frames(list_of_frames, outfile, iheader=0):
    """Compute averaged frame from list of frames.
    
    Parameters
    ----------
    list_of_frames : list of strings
        List containing the list of invidual FITS files to be avaraged.
    outfile : string
        Outpuf file name
    iheader : int
        Image number (in list_of_frames) from which the header will
        be copied. If 0, no header is added to the output file.
        
    """

    nframes = len(list_of_frames)

    main_header = None

    for i in range(nframes):
        infile = list_of_frames[i]
        print("\n*Working with file:\n", infile)
        hdulist = fits.open(infile)
        if i == iheader:
            main_header = hdulist[0].header
        image2d = hdulist[0].data.astype(np.float)
        hdulist.close()
        if i == 0:
            result = np.copy(image2d)
        else:
            result += image2d

    result /= nframes
    if main_header is None:
        hdu = fits.PrimaryHDU(result)
    else:
        hdu = fits.PrimaryHDU(result, main_header)  # add main_header
    hdu.writeto(outfile, clobber=True)
    print("\n*Generating output file:\n" + outfile)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog='reduce')
    parser.add_argument("input_darks", 
                        help="txt file with list of initial darks")
    parser.add_argument("input_objects", 
                        help="txt file with list of initial images")

    parser.add_argument("nramps",
                        help="number of ramps",type=int)
    args = parser.parse_args()
    
    list_darks = np.genfromtxt(args.input_darks, dtype=[('darksname', '|S100')])
    darksname = list_darks['darksname']
    
    list_objects = np.genfromtxt(args.input_objects, dtype=[('objectsname', '|S100')])
    objectsname = list_objects['objectsname']

    nramps = int(args.nramps)
    
    
    execute_darks = True
    # compute averaged dark
    if execute_darks:
        list_darks_coadd = []
        for infile in darksname:
            outfile = infile[:-5] + "_coadd.fits"
            list_darks_coadd.append(outfile)
            read_raw_circe(infile, nramps, outfile)
        average_frames(list_darks_coadd, "dark.fits", iheader=0)

    execute_objects = True
    
    if execute_objects:
        
        list_object_coadd = []
        for infile in objectsname:
            outfile = infile[:-5] + "_coadd.fits"
            list_object_coadd.append(outfile)
            read_raw_circe(infile, nramps, outfile, debug=1)




