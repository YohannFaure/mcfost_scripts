#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from mcfostModule import *


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'arguments')
    parser.add_argument("loc1", help = "file with the beam properties in it's header.",type = str)
    parser.add_argument("loc2", help = "MCFOST File.",type = str)
    parser.add_argument("output", help = "Output File",type = str)
    parser.add_argument("--fastplot", help='If you want to make a fast plot of the result.',action='store_true')
    args = parser.parse_args()
    ##### Get args
    loc1,loc2=args.loc1,args.loc2
    output=args.output
    ##### Make beam
    beam=MakeBeamMesh(loc1,loc2)
    ##### Open The two images
    _,Im1=openimage(loc1)
    header,Im2=openimage(loc2)
    bmin,bmaj,bpa=ExtractBeam(loc1)

    # A full convolution is slow, I need a smaller grid.
    l=len(beam)
    maxbeam=np.max(beam[l//2])
    i_lim=0
    while beam[l//2,i_lim]<maxbeam/100:
            i_lim+=1
    nconv=l//2-i_lim
    #print(nconv)
    ##### Convolution
    beamsmall = beam[l//2-nconv:l//2+nconv,l//2-nconv:l//2+nconv]
    from scipy.signal import fftconvolve
    ConvolvedImage=fftconvolve(Im2,beamsmall)

    ##### If you want a plot
    if args.fastplot:
        fig,axl=plt.subplots(2,2)
        ax=axl.flatten()
        ax[0].imshow(ConvolvedImage)
        ax[1].imshow(Im2)
        ax[2].imshow(Im1)
        ax[3].imshow(beam)
        plt.show()

    ##### Saves the data
    header.set('BPA',bpa)
    header.set('BMIN',bmin)
    header.set('BMAJ',bmaj)
    header.set('IMAGE_BY','Yohann Faure')
    header.set('BEAMTYPE', 'Gaussian')

    CreateFitsFile(output,header,ConvolvedImage)

