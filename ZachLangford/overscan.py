#!/usr/bin/python

# Version 2.1 -- Elinor Gates, 2016 Aug 25 - changed math to BITPIX -32
# Version 2.0 -- Elinor Gates, 2016 Aug 11 - improved file reading
# Version 1.0 -- Elinor Gates, 2015 Nov 24


from astropy.io import fits,ascii
import numpy as np
import sys, getopt

def main(argv):
    inputfilelist = ''
    outputfilelist = ''
    fit = 'no'
    try:
        opts, args = getopt.getopt(argv,"hfi:o:",["ifilelist=","ofilelist="])
    except getopt.GetoptError:
      print 'overscanLickObs.py -f -i <inputfilelist> -o <outputfilelist>'
      print '-f indicates do a Legendre fit to overscan'
      sys.exit(2)
    for opt, arg in opts:
      if opt == '-h':
         print 'overscanLickObs.py -f -i <inputfilelist> -o <outputfilelist>'
         print '-f indicates do a Legendre fit to overscan'
         sys.exit(2)
      elif opt in ("-i", "--ifilelist"):
         inputfilelist = arg
      elif opt in ("-o", "--ofilelist"):
         outputfilelist = arg
      elif opt == '-f':
          fit = 'yes' 
    print 'Input filelist is ', inputfilelist
    print 'Output filelist is ', outputfilelist
    print 'Fit is ', fit

    # open input and output filelists
    ifilelist = [line.rstrip('\n') for line in open(inputfilelist)]
    ofilelist = [line.rstrip('\n') for line in open(outputfilelist)]
    
    # how many files
    numifiles = len(ifilelist)
    numofiles = len(ofilelist)
    if numifiles != numofiles:
        sys.exit('Input and output file lists have different numbers of files. Exiting.')

    # for each file in ifilelist, read in file, figure out overscan and data regions, fit 
    # overscan with desired function (if any), and subtract from data.  
    # Write data to ofilelist value.  

    for i in range(0,numifiles):
        ifile=ifilelist[i]
        ofile=ofilelist[i]
        data, header = fits.getdata(ifile,header=True)

        # change data to float
        data=data.astype('float32')

        # read necessary keywords from fits header
        xsize = header['NAXIS1']
        ysize = header['NAXIS2']
        xorig = header['CRVAL1U']
        yorig = header['CRVAL2U']
        cdelt1 = header['CDELT1U']
        cdelt2 = header['CDELT2U']
        rover = header['ROVER']
        cover = header['COVER']
        inxsize = header['DNAXIS1']
        inysize = header['DNAXIS2']
        ampsx = header['AMPSCOL']
        ampsy = header['AMPSROW']

        # determine number and sizes of overscan and data regions
        namps = ampsx*ampsy
        if rover > 0:
            over=rover
            sys.exit('Program does not yet deal with row overscans. Exiting.')
        else:
            over = cover
        if over == 0:
            sys.exit('No overscan region specified in FITS header. Exiting.')

        # single amplifier mode
        if namps == 1:
            biassec = data[:,xsize-cover:xsize]
            datasec = data[0:,0:xsize-cover]

            # median overscan section
            bias=np.median(biassec, axis=1) 

            # legendre fit
            if fit == 'yes':
                # fit
                lfit = np.polynomial.legendre.legfit(range(0,len(bias)),bias,3)
                bias = np.polynomial.legendre.legval(range(0,len(bias)),lfit)

            # subtract overscan
            datanew = datasec
            for i in range(datasec.shape[1]):
                datanew[:,i] = datasec[:,i]-bias

        # two amplifier mode
        if namps == 2:
            biasseca = data[:,xsize-cover*2:xsize-cover]
            biassecb = data[:,xsize-cover:xsize]

            # median overscan sections
            biasa=np.median(biasseca,axis=1)
            biasb=np.median(biassecb,axis=1)

            # legendre fit
            if fit == 'yes':
                lfita = np.polynomial.legendre.legfit(range(0,len(biasa)),biasa,3)
                lfitb = np.polynomial.legendre.legfit(range(0,len(biasb)),biasb,3)
                biasa = np.polynomial.legendre.legval(range(0,len(biasa)),lfita)
                biasb = np.polynomial.legendre.legval(range(0,len(biasb)),lfitb)

            # extract data regions

            #determine size of binned data region
            hsize=abs(inxsize/cdelt1)
 
            # calculate x origin of readout in binned units if cdelt1 negative or positive
            if cdelt1 < 0:
                xorig=(xorig-(xsize-2*cover)*abs(cdelt1))/abs(cdelt1)
            else:
                xorig=xorig/cdelt1
            x0=xorig+xsize-1-cover*2 # need to test is need -1 because starting counting at 0

            # determine which columns are on which amplifier and subtract proper overscan region

            if x0 < hsize/2: # all data on left amplifier
                datanew=data[:,0:xsize-cover*2]
                m=datanew.shape[1]
                for i in range(0,m):
                    datanew[:,i]=datanew[:,i]-biasa

            if xorig >= hsize/2: # all data on right amplifier
                datanew=data[:,0:xsize-cover*2]
                m=datanew.shape[1]
                for i in range(0,m):
                    datanew[:,i]=datanew[:,i]-biasb

            if xorig < hsize/2 and x0 > hsize/2:
                x1=hsize/2-xorig
                dataa=data[:,0:x1]
                datab=data[:,x1:-cover*2]
                ma=dataa.shape[1]
                mb=datab.shape[1]
                for i in range(0,ma):
                    dataa[:,i]=dataa[:,i]-biasa
                for i in range(0,mb):
                    datab[:,i]=datab[:,i]-biasb
                # merge dataa and datab into single image
                datanew=np.hstack([dataa,datab])

        if namps > 2: 
            sys.exit('Program does not yet deal with more than two overscan regions. Exiting.')

        # add info to header
        header['HISTORY'] = 'Overscan subtracted'

        # write new fits file
        fits.writeto(ofile,datanew,header,clobber=True)

if __name__ == "__main__":
   main(sys.argv[1:])




