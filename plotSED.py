#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from mcfostModule import *



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'arguments')
    parser.add_argument("sed_folder", help = "Loaction of the SED.",type = str)
    parser.add_argument("param_location", help = "MCFOST para File.",type = str)
    parser.add_argument("--XPSED", help = "Experimental SED",action='store_true')
    parser.add_argument("--Av", help="extinction", type=float, default=None)
    #parser.add_argument("--fastplot", help='If you want to make a fast plot of the result.',action='store_true')
    args = parser.parse_args()
    ##### Define the working folder
    folder=args.sed_folder+'/'
    filename='sed_rt.fits.gz'
    parameters=args.param_location

    ##### Extract the data
    datalist,param=opendata(filename,parameters,folder)
    sed_profile = np.squeeze(datalist[0].data)
    wl=WaveLengthExtract(param,folder)

    ##### extinct
    if args.Av:
        sed_profile = extinct(args.Av,wl,sed_profile)

    ##### Plot
    plotSED(wl,sed_profile,label='MCFOST data.',color='k')

    if args.XPSED:
        ##### Same with XP
        lsed=['SEDPoints/SED_J1615_dict.npy','SEDPoints/SED_KOOISTRA_dict.npy','SEDPoints/SED_MAREL_dict.npy']#,'SEDPoints/SED_VOSA_dict.npy']
        captions=[r'Our Alma observation',r'Kooistra et al. 2016',r'Van der Marel et al. 2016']#,r'VOSA Database']
        color=['g','b','r']#,'y']
        #wl,sed_profile,error=SEDfromVOSA(datfile)
        #sed_profile=ergscm2aTOwm2(sed_profile,wl)
        wl,sed,err,origin=ExperimentalSED(lsed)
        for i in range(len(wl)):
            plotSED(wl[i],sed[i],err[i],label=captions[i],color=color[i])
    plt.legend()
    plt.savefig('SED.pdf')
