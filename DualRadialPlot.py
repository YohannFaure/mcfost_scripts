#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot a superposition of two radial plots

python3 DualRadialPlot.py result.fits J1615_edit.fits 45.9835388 326.2477693  6.341708607533389 3. --CenterAndWidth1 '(1555,1555,600)' --CenterAndWidth2 '(1500,1500,600)' --label1 'MCFOST Simulation' --label2 'Alma Data'
"""

from mcfostModule import *



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'arguments')
    parser.add_argument("im1", help = "First Image.",type = str)
    parser.add_argument("im2", help = "Second image.",type = str)
    parser.add_argument("inc", help = "inclination",type=float)
    parser.add_argument("pa", help = "position angle.",type = float)
    parser.add_argument("par", help = "paralax.",type = float)
    parser.add_argument("binwidth", help = "Binning width.",type = float)
    parser.add_argument("--CenterAndWidth1", help = "First Image center and width '(x,y,w)'.",type = str,default=None)
    parser.add_argument("--CenterAndWidth2", help = "Second Image center and width '(x,y,w)'.",type = str,default=None)
    parser.add_argument("--label1", help = "First image label",type = str,default=None)
    parser.add_argument("--label2", help = "Second image label",type = str,default=None)
#    parser.add_argument("--Av", help="extinction", type=float, default=None)
    #parser.add_argument("--fastplot", help='If you want to make a fast plot of the result.',action='store_true')
    args = parser.parse_args()
    r1,m1,s1=NormalizedRadialProfile(args.im1,args.inc,args.pa,args.binwidth,par=args.par,CenterAndWidth=args.CenterAndWidth1,save=None)
    r2,m2,s2=NormalizedRadialProfile(args.im2,args.inc,args.pa,args.binwidth,par=args.par,CenterAndWidth=args.CenterAndWidth2,save=None)
    fig=plt.figure()
    fig.set_size_inches(5,4)
    PlotRadialProfile(r1,m1,s1,unit='normalized',color='r',label=args.label1)
    PlotRadialProfile(r2,m2,s2,unit='normalized',color='b',label=args.label2)
    plt.legend()
    plt.show()
    #input('Press enter to close.')
