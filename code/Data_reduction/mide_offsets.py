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
import argparse
import matplotlib.pyplot as plt
import numpy as np
import re
from skimage.feature import register_translation

from pause_debugplot import pause_debugplot


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


def cross_correlate(ref_image, shifted_image, debugplot=0):
    """Cross-correlation to identify the relative shift.

    Parameters
    ----------
    ref_image : 2d numpy array (float)
        Reference image.
    shifted_image : 2d numpy array (float)
        Shifted image.
    debugplot : int
        Determines whether intermediate computations and/or plots
        are displayed:
        00 : no debug, no plots
        01 : no debug, plots without pauses
        02 : no debug, plots with pauses
        10 : debug, no plots
        11 : debug, plots without pauses
        12 : debug, plots with pauses

    Returns
    -------
    shifts : ndarray
        Shift vector (in pixels) required to register ``shifted_image``
         with ``reference``. Axis ordering is consistent with numpy
         (e.g. Z, Y, X).
    error : float
        Translation invariant normalized RMS error between
        ``reference`` and ``shifted_image``.
    phasediff : float
        Global phase difference between the two images (should be
        zero if images are non-negative).

    """

    shifts, error, diffphase = \
        register_translation(ref_image, shifted_image)

    if debugplot >= 10:
        print(">>> Offset (y, x):", shifts)
        print(">>> Error........:", error)
        print(">>> Diffphase....:", diffphase)

    if debugplot % 10 != 0:
        fig = plt.figure(figsize=(12, 5))
        ax1 = plt.subplot(1, 3, 1, adjustable='box-forced')
        ax2 = plt.subplot(1, 3, 2, sharex=ax1, sharey=ax1,
                          adjustable='box-forced')
        ax3 = plt.subplot(1, 3, 3, sharex=ax1, sharey=ax1,
                          adjustable='box-forced')

        ax1.imshow(ref_image, interpolation='nearest', origin='low')
        # ax1.set_axis_off()
        ax1.set_title('Reference image')

        ax2.imshow(shifted_image, interpolation='nearest', origin='low')
        # ax2.set_axis_off()
        ax2.set_title('Offset image')

        # View the output of a cross-correlation to show what the
        # algorithm is doing behind the scenes
        image_product = np.fft.fft2(ref_image) * \
                        np.fft.fft2(shifted_image).conj()
        cc_image = np.fft.fftshift(np.fft.ifft2(image_product))
        ax3.imshow(cc_image.real, interpolation='nearest', origin='low')
        # ax3.set_axis_off()
        ax3.set_title("Cross-correlation")

        plt.show(block=False)
        plt.pause(0.001)
        pause_debugplot(debugplot)
        plt.close()

    return shifts, error, diffphase


if __name__ == "__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument("input_list", 
                        help="txt file with list of initial images")
    parser.add_argument("input_offsets", 
                        help="txt file (from xnirspec) with list of " + \
                        "initial offsets")
    parser.add_argument("x0",
                        help="x0 coordinate corresponding to object in image")
    parser.add_argument("y0",
                        help="y0 coordinate corresponding to object in image")
    parser.add_argument("delta_x",
                        help="x-width of the image subregion to be employed")
    parser.add_argument("delta_y",
                        help="y-width of the image subregion to be employed")
    parser.add_argument("--dont_subtract_consecutive",
                        help="avoid subtraction of  consecutive images",
                        action="store_true")
    parser.add_argument("--dark_name",
                        help="file name of the dark file (default: dark.fits)",
                        default="dark.fits")
    parser.add_argument("--debugplot",
                        help="integer indicating plotting/debugging" + \
                        " (default=0)", type=int,
                        default=0)
    args = parser.parse_args()

    x0 = int(args.x0)
    y0 = int(args.y0)
    delta_x = int(args.delta_x)
    delta_y = int(args.delta_y)

    list_files = np.genfromtxt(
        args.input_list, dtype=[('filename', '|S100')])
    list_offsets = np.genfromtxt(
        args.input_offsets, 
        dtype=[('iflag', '<i8'), ('filename', '|S100'),
               ('offx', '<i8'), ('offy', '<i8')]
    )
    filename = list_offsets['filename']
    #filename = list_files['filename']
    offx = list_offsets['offx']
    offy = list_offsets['offy']
    nfiles = filename.size

    for val in zip(filename, offx, offy):
        print(val)

    list_regions = []
    for k in range(nfiles):
        if k < nfiles - 1:
            image1, dum, dum = readfits(filename[k])
            image2, dum, dum = readfits(filename[k + 1])
        else:
            image1, dum, dum = readfits(filename[k])
            image2, dum, dum = readfits(filename[k - 1])
            
        if args.dont_subtract_consecutive:
            image = image1
        else:
            image = image1 - image2


        i1 = int(y0 - delta_y/2 - offy[k])
        i2 = int(y0 + delta_y/2 - offy[k])
        j1 = int(x0 - delta_x/2 - offx[k])
        j2 = int(x0 + delta_x/2 - offx[k])
        sublist_regions = []
        npixstep=2  # step size (in pixels)
        number_of_semisteps = 2  # at both sides of the initial position
        for deltai in np.arange(-number_of_semisteps, number_of_semisteps+1, 1, dtype=int)*npixstep:
            for deltaj in np.arange(-number_of_semisteps, number_of_semisteps+1, 1, dtype=int)*npixstep:
                region = image[(i1+deltai):(i2+deltai), (j1+deltaj):(j2+deltaj)]
                sublist_regions.append(region)
        list_regions.append(sublist_regions)

    sublist_reference = list_regions[0]
    for k in range(1,nfiles):
        sublist_offset_region = list_regions[k]
        if args.debugplot >= 10:
            print("* Reference image:", filename[0])
            print("* Offset image...:", filename[k])
        l = 0
        offx_tmp = np.zeros((2*number_of_semisteps+1)*(2*number_of_semisteps+1))
        offy_tmp = np.zeros((2*number_of_semisteps+1)*(2*number_of_semisteps+1))
        for deltai in np.arange(-number_of_semisteps, number_of_semisteps+1 ,1, dtype=int)*npixstep:
            for deltaj in np.arange(-number_of_semisteps, number_of_semisteps+1, 1, dtype=int)*npixstep:
                tmp_reference = sublist_reference[l]
                tmp_offset_region = sublist_offset_region[l]
                shifts, error, diffphase = \
                    cross_correlate(tmp_reference, tmp_offset_region, debugplot=args.debugplot)
                offx_tmp[l] = shifts[1]
                offy_tmp[l] = shifts[0]
                print('>>>', l, shifts[0], shifts[1])
                l += 1
        print(offx_tmp)
        print(offy_tmp)
        print('>>> median offx, std offx: ', np.median(offx_tmp), np.std(offx_tmp))
        print('>>> median offy, std offy: ', np.median(offy_tmp), np.std(offy_tmp))
        offy[k] += np.median(offy_tmp)
        offx[k] += np.median(offx_tmp)
        print(k, offx[k], offy[k])


    print("--- ooo ---")
    for val in zip(filename, offx, offy):
        print(val)

    list_offsets['filename'] = list_files['filename']
    output_txt = re.sub('.txt$', '_refined.txt', args.input_offsets)
    np.savetxt(output_txt, list_offsets,
               fmt="%1d %s " + args.dark_name + " %6d %6d")
