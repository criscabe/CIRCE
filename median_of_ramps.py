from __future__ import division
from __future__ import print_function
from astropy.io import fits
import argparse
import numpy as np
import re
import sys


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
    
    # Only to ckeck the images dimensions 
    NAXIS1 = main_header.get('NAXIS1')
    NAXIS2 = main_header.get('NAXIS2')
    # print(NAXIS1) 
    # print(NAXIS2)
    
    hdulist.close()

    return image2d, main_header, header
   

def trimming_images(image2, offsetx, offsety):
      
      
      """Trimming of images.

    Parameters
    ----------
    
    image2 : 2d numpy array (float)
        Successive ramps.(Not the first)
        
        
    offsetx: int 
        Offset in x between image2 and image1
        
    offysety: int 
        Offset in y between image2 and image1

    Returns
    -------
    IMP2 : 2d numpy array, floats
        Data array corresponding to trimmed image

    """ 
    
      i1 = j1 = 0 
      i2 = image2.shape[0]
      j2 = image2.shape[1]

      IMP2= np.zeros((i2,j2)) 
      
      if offsetx >= 0 and offsety > 0 :
          IMP =  image2[(i1+offsety):i2, j1:(j2-offsetx)] 
          IMP2[:(i2-offsety), (j1+offsetx):] = IMP

      elif offsetx <= 0 and offsety < 0 :
          IMP=  image2[i1:(i2+offsety), (j1-offsetx):j2]
          IMP2[(i1-offsety):,: (j2+offsetx)]=IMP

      elif offsetx < 0 and offsety >= 0 :
          IMP=  image2[(i1+offsety):i2, (j1-offsetx):j2]
          IMP2[:(i2-offsety), : (j2+offsetx)]=IMP
      
      elif offsetx > 0 and offsety <= 0 :
          IMP=  image2[i1:(i2+offsety), j1:(j2-offsetx)]
          IMP2[(i1-offsety):, (j1+ offsetx):]=IMP

         
      else :  # both offsets are zero
          IMP2 = np.copy(image2)
         
    
      return IMP2
    
if __name__ == "__main__":
    

    # define expected command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_list", 
                 help='"txt file with list of initial file names. Initial ramps"')

    parser.add_argument("input_offsets", 
                        help='"txt file  with list of offsets between ramps"')
    
    
    parser.add_argument("--cube_with_ramps",
                        help="generate cube with ramps",
                        action="store_true")
    parser.add_argument("--trimmed_images",
                        help="generate the trimmed ramps",
                        action="store_true")
    
    args = parser.parse_args()

    
    # read list with file names
    
    list_files = np.genfromtxt(args.input_list, dtype=[('filename_initial', '|S100')])
    filename_initial = list_files['filename_initial']

    nfiles_initial = filename_initial.size
    # nfiles_initial = number of ramps

    list_offsets = np.genfromtxt(args.input_offsets, 
        dtype=[('iflag', '<i8'), ('filename', '|S100'),('dark', '|S100'),
               ('offx', '<i8'), ('offy', '<i8')] )
    
    
    filename = list_offsets['filename']

    offx = list_offsets['offx']
    offy = list_offsets['offy']

    # Print the offsets between ramps
    for val in zip(filename, offx, offy):
        print(val)
        

    if offx[0] != 0 or offy[0] != 0:
        sys.exit("ERROR!!: Reference image has offsets different to zero")
    
    naxis1 = 2048
    naxis2 = 2048
        
    image3d = np.zeros((nfiles_initial , naxis2, naxis1))  

    for k in range(nfiles_initial):
        image, main_header, dum = readfits(filename_initial[k])
        if image.shape != (naxis2, naxis1):
            sys.exit("ERROR!!: unexpected image dimensions")
            
        if k == 0:
            image3d[0,:,:] = image[:,:]
        else:        
            T_image = trimming_images(image, offx[k], -offy[k] )
            image3d[k,:,:] = T_image[:,:]

    median = np.median(image3d, axis=0)
    
    #  Generating the median image of the ramps 
    
    outfile2 = re.sub('ramp1.fits$', 'median.fits', filename_initial[0])
    print("Generating the median of the ramps:\n" + outfile2)  
    hdu = fits.PrimaryHDU(median, main_header)
    hdu.writeto(outfile2, clobber=True)

    #  Generating cube with the ramps (corrected the offsets)
    #  In principle we dont need generate the cube.fits
    
    if args.cube_with_ramps:
        outfile = re.sub('ramp1.fits$', 'cube_with_ramps.fits', filename_initial[0])
        print("Generating datacube with all the ramps:\n" + outfile)
        hdu = fits.PrimaryHDU(image3d, main_header)
        hdu.writeto(outfile, clobber=True)
    
    #  Generating the trimmed images of the ramps.
    #  In principle we dont need generate the trimmed images
    
    if args.trimmed_images:
        list_cube_split = []  
        for k in range(nfiles_initial):
            list_cube_split.append(image3d[k,:,:])
            outfile1 = re.sub('.fits$', '_trimmed.fits', filename_initial[k])
            print("Generating trimmed images:\n" + outfile1)
            hdu = fits.PrimaryHDU(list_cube_split[k], main_header)
            hdu.writeto(outfile1, clobber=True)
        