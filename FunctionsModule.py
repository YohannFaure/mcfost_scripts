#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Basinc data fitting functions for astrophysics.
Python3.7

Yohann Faure
ENS de Lyon
yohann.faure@ens-lyon.fr
'''

import sys
import numpy as np
import matplotlib.pyplot as plt
#sys.path.append('/home/nicolas/Documents/Programas/astronomia/trabajos/2019/CSCha')
from astropy.io import fits
from matplotlib.colors import LogNorm
from scipy.interpolate import interp1d
from astropy.visualization import (AsinhStretch, LogStretch, LinearStretch, ImageNormalize)
from scipy import ndimage
import scipy.optimize as opt
from time import time
import multiprocessing as mp
import matplotlib as mpl
from scipy.signal import savgol_filter


radtodeg = 180./np.pi
degtorad=np.pi/180.

##############################################################################
#                                   FUNCTIONS                                #
##############################################################################

##### Geometry

def deproject(x, y, inc, pa):
    '''
    Takes positions x and y, and deprojects them usins the inclination
    and position angle given, in radians.
    '''
    cosinc=np.cos(inc)
    sinpa=np.sin(pa)
    cospa=np.cos(pa)
    xp = x * (cospa/cosinc) + y * (sinpa/cosinc)
    yp = y * cospa - x * sinpa
    return xp, yp

def project(x, y, inc, pa):
    '''
    Takes positions x and y, and projects them usins the inclination
    and position angle given, in radians.
    '''
    cosinc=np.cos(inc)
    sinpa=np.sin(pa)
    cospa=np.cos(pa)
    xpi = x * cosinc
    xp = xpi * cospa - y * sinpa
    yp = y * cospa + xpi * sinpa
    return xp, yp

def deprojectedcoordinates(x,y,inc,pa):
    xx,yy=deproject(x,y,inc,pa)
    rr=np.sqrt(xx**2+yy**2)
    angles=np.arctan2(xx,yy)
    return(xx,yy,rr,angles)

def deprojectedcoordinatespairs(pair,inc,pa):
    x,y=pair
    xx,yy=deproject(x,y,inc,pa)
    rr=np.sqrt(xx**2+yy**2)
    angles=np.arctan2(xx,yy)
    return(xx,yy,rr,angles)

def angular_difference(i, j):
    '''
    Calculates the difference between two angles (deg).
    '''
    diff = ((i - j + 180.) % 360.) - 180. #Parkour
    return diff

def gaussiana(r, A, r0, sigma):
    '''
    Calculates the gaussian value of the position r, given that the gaussian
    is centered at gr0 and has a 1 sigma width of gsig. The gaussian amplitude
    is given by gf0.
    '''
    profile = A * np.exp(-0.5 * ( (r - r0) / sigma )**2)
    return profile

def gaussianfit(param, mesh):
    """
    create the gaussian profile associated with th eimage, with parameters
    x0, y0 : center
    inc, pa: inclination and position angle
    a, sigma = amplitude and sigma of the gaussian
    mesh = meshgrid for the fitting. Computing it out of the functions allows a faster minimization.
    """
    # Params
    x0, y0, inc, pa, a, sigma = param
    # Calculate pixel grid
    x,y = mesh
    x,y=x-x0,y-y0
    # Deproject distances
    dep_x, dep_y = deproject(x, y, inc, pa)
    dep_r = np.sqrt(dep_x**2 + dep_y**2)
    return(gaussiana(dep_r, a, 0., sigma))

def radialgaussian(param,mesh):
    A,r0,sigma,inc,pa = param
    x,y=mesh
    xx, yy = deproject(x, y, inc, pa)
    rr = np.sqrt(xx**2 + yy**2)
    return(gaussiana(rr, A, r0, sigma))

##### Imaging

def meshmaker(llx,lly,xc,yc):
    x = np.arange(llx)-xc
    y = np.arange(lly)-yc
    return(np.meshgrid(x, y))

def openimage(location):
    """
    opens the image and the header
    """
    hdulist = fits.open(location)
    header = hdulist[0].header
    image = np.squeeze(hdulist[0].data)
    return(header, image)

def resizedimage(image,xc,yc,xsize,ysize=None):
    """
    resizes an image
    xc,yc=center coordinates
    xsize,ysize= new size of the image
    """
    if ysize == None: ysize=xsize
    x_in, x_fi, y_in, y_fi = xc - xsize, xc + xsize, yc - ysize, yc + ysize
    image = image[x_in:x_fi, y_in:y_fi]
    return(image)

def extractheaderinfo(header):
    """
    extracts some info, you should prefer the direct way: x=header['...']
    """
    ra  = 3600. * header['CDELT1'] * (np.arange(header['NAXIS1']) - (header['CRPIX1']-1))
    dec = 3600. * header['CDELT2'] * (np.arange(header['NAXIS2']) - (header['CRPIX2']-1))
    ext = [np.max(ra), np.min(ra), np.min(dec), np.max(dec)]
    return(ra,dec,ext)

def blurimage(img,sigma=30):
    """Gaussian-blurs an image"""
    blur = ndimage.gaussian_filter(img, sigma=sigma)
    return(blur)

def GlobalVars(header):
    pixtosec=np.abs(3600*header['CDELT1'])
    par=header['PAR']
    pixtoau=1000*pixtosec/(par)
    arcsectoau=1000/par
    return(pixtosec,pixtoau,arcsectoau)

def deprojectimage(location,inc,pa,resizeparams):
    """
    Makes an image of a tilted disk as if facing us.
    resizeparams coresponds to the parameters in input of the "resize" function. It is higly recomended to crop the image as this function is a bit slow.
    """
    ##### Open initial image and get info
    header,image=openimage(location)
    lx,ly,w=resizeparams
    par=header["PAR"]
    pixtosec=np.abs(3600*header['CDELT1'])
    pixtoau=1000*pixtosec/(par)
    arcsectoau=1000/par
    image=resizedimage(image,lx,ly,w)
    ##### Make mesh
    llx,lly=image.shape
    x, y = meshmaker(llx,lly,llx/2,lly/2)
    xx,yy=deproject(x,y,inc*degtorad,pa*degtorad)
    rp=np.sqrt(xx**2+yy**2)
    ##### Other way
    #from mpl_toolkits.mplot3d import Axes3D
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.plot_surface(xx, yy, image,cmap='inferno',interpolation=nearest)
    #plt.show()
    ##### Plot in itself
    fig = plt.figure()
    fig.set_size_inches(5,5)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_facecolor('grey')
    fig.suptitle('Our protoplanetary disk front-viewed',size=12)
    ##### see documentation for more info. This is the best I found.
    ax.contourf(xx*pixtoau,yy*pixtoau,image,cmap='inferno',levels=1000,interpolation='nearest')
    ax.axis('equal')
    ax.set_xlabel(r"$\vec{u_x}$ (au)")
    ax.set_ylabel(r"$\vec{u_y}$ (au)")
    ##### Real units and beam
    plt.xlim(np.max(xx)*pixtoau,np.min(xx)*pixtoau)
    beam = mpl.patches.Ellipse(xy=(2.50*arcsectoau,-2.50*arcsectoau), width=header['BMIN']*3600*arcsectoau, height=header['BMAJ']*3600*arcsectoau , color='white', fill=True, angle=header['BPA'])
    ax.add_artist(beam)
    plt.savefig('frontview{:.0f}.png'.format(time()),dpi=600)
    return(None)

def deprojectimageBIS(header,image,inc,pa):
    """
    Makes an image of a tilted disk as if facing us.
    resizeparams coresponds to the parameters in input of the "resize" function. It is higly recomended to crop the image as this function is a bit slow.
    """
    ##### Open initial image and get info
    pixtosec=np.abs(3600*header['CDELT1'])
    ##### Make mesh
    llx,lly=image.shape
    x, y = meshmaker(llx,lly,llx/2,lly/2)
    xx,yy=deproject(x,y,inc*degtorad,pa*degtorad)
    rp=np.sqrt(xx**2+yy**2)
    ##### Other way
    #from mpl_toolkits.mplot3d import Axes3D
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.plot_surface(xx, yy, image,cmap='inferno',interpolation=nearest)
    #plt.show()
    ##### Plot in itself
    fig = plt.figure()
    fig.set_size_inches(5,5)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_facecolor('grey')
    fig.suptitle('Our protoplanetary disk front-viewed',size=12)
    ##### see documentation for more info. This is the best I found.
    ax.contourf(xx*pixtosec,yy*pixtosec,image,cmap='inferno',levels=1000,interpolation='nearest')
    ax.axis('equal')
    ax.set_xlabel(r"$\vec{u_x}$ (arcsec)")
    ax.set_ylabel(r"$\vec{u_y}$ (arcsec)")
    ##### Real units and beam
    plt.xlim(np.max(xx)*pixtosec,np.min(xx)*pixtosec)
    beam = mpl.patches.Ellipse(xy=(0.9*np.max(xx)*pixtosec,-0.9*np.max(xx)*pixtosec), width=header['BMIN']*3600, height=header['BMAJ']*3600 , color='white', fill=True, angle=header['BPA'])
    ax.add_artist(beam)
    plt.savefig('frontview{:.0f}.png'.format(time()),dpi=600)
    return(None)

##### Misc

def multithreadmap(f,X,ncores=8):
	"""
	multithreading map of a function, default on 8 cpu cores.
    might be usefull someday.
	"""
	p=mp.Pool(ncores)
	Xout = p.map(f,X)
	p.terminate()
	return(Xout)

def slidemean(a,nslide):
    """sliding mean of an array"""
    return(np.convolve(a, np.ones((nslide,))/nslide, mode='valid'))

def polynomean(a,nslide,order):
    """Polynomial sliding mean"""
    return(savgol_filter(a, nslide, order))

def quadcost(img,model):
    return(np.sum((img-model)**2))

def showi(image):
    plt.imshow(image)
    plt.show()

##### Profiles

def radialbin(image, radii, binwidth):
    """
    Makes a radial binning of an image with given radii, with a given bin width
    """
    radiusbinned=(radii//binwidth).astype(int) #chosing the bin
    n=np.max(radiusbinned) # to set rmax and rmin
    I_of_r=[[] for i in range(n+1)]
    l1,l2=image.shape
    for x in range(l1):
        for y in range(l2):
            I_of_r[radiusbinned[x,y]].append(image[x,y])
    means=[]
    stds=[]
    for k in I_of_r:
        means.append(np.mean(k))
        stds.append(np.std(k))
    return(means,stds)

def angularprofile(image,rr,angles,radius,rwidth):
    """
    returns an angular profile for an image with computed radii rr and angles, at fixed radius, with a tolerence of rwidth
    image = the image array
    rr = array of every radius in the image (same shape as image)
    angles = array of angles (same shape too)
    radius = radius at which the profile is made
    rwidth = tolerence
    """
    coords=[]
    l1,l2=image.shape
    truth = (rr<radius+rwidth) * (rr> radius-rwidth)
    for i in range(l1):
        for j in range(l2):
            if truth[i,j]:
                coords.append([angles[i,j],image[i,j]])
                image[i,j]=None
    return(np.array(coords),image)


##### Higly custom, probably not what you are looking for
def modelplot2(model,image,cmap,header,size=(4.77,4.3),save='QuadplotSecondEdition{}.png'.format('test')):#time())):
    """
    Plots a quadri-images graph, need to be adapted to your own plot, as it is higgly custom...
    """
    ##### settings
#    if save!='False':
#        mpl.use('pgf')
#        pgf_with_rc_fonts = {
#            "font.family": "serif",
#            "font.serif": [],                   # use latex default serif font
#            "font.sans-serif": [], # use a specific sans-serif font
#            "font.size": 8}
#        plt.style.use('seaborn-whitegrid')
#        mpl.rcParams.update(pgf_with_rc_fonts)
    maximage=np.max(image)
    image,model=image/maximage,model/maximage
    fig,axarr=plt.subplots(2,2)
    plt.suptitle('Fitting of a brightness profile\n(normalization: {:5.5f} {})'.format(maximage,header['BUNIT']),fontsize=11)
    fig.set_size_inches(size[0],size[1])
    ##### scale for pixel--arcsec conversion
    par=header["PAR"]
    pixtosec=np.abs(3600*header['CDELT1'])
    pixtoau=1000*pixtosec/(par)
    arcsectoau=1000/par
    scalex=image.shape[0]*pixtoau/2
    scaley=image.shape[1]*pixtoau/2
    ##### setting and plotting each subplots
    ### Titles
    axarr[0,0].title.set_text(r'Fitting model')
    axarr[0,1].title.set_text(r'10$\times$Residual')
    axarr[1,0].title.set_text(r'Original image')
    axarr[1,1].title.set_text(r'Residual\\absolute value')
    ### Plots
    axarr[0,0].imshow(model, vmin=0, vmax=1, interpolation='nearest',origin=0,cmap=cmap,extent=[scalex,-scalex,-scaley, scaley])
    axarr[0,1].imshow(10*np.abs(model-image), vmin=0, vmax=1, interpolation='nearest',origin=0,cmap=cmap,extent=[scalex,-scalex,-scaley, scaley])
    axarr[1,0].imshow(image, vmin=0, vmax=1, interpolation='nearest',origin=0,cmap=cmap,extent=[scalex,-scalex,-scaley, scaley])
    ### resolution bean on the 3rd subplot
    beam = mpl.patches.Ellipse(xy=(1.75*arcsectoau,-1.75*arcsectoau), width=header['BMIN']*3600*arcsectoau, height=header['BMAJ']*3600*arcsectoau , color='white', fill=True, angle=header['BPA'])
    axarr[1,0].add_artist(beam)
    ### legend
    axarr[1,0].set_xlabel('East-West (au)')
    axarr[1,0].set_ylabel('South-North (au)')
#    axarr[1,1].set_xlabel('Inclination: {inc:3.2f}. Position angle: {pa:3.2f}. (deg)'.format(inc=params[2]*radtodeg,pa=params[3]*radtodeg))
    ### last subplot
    im=axarr[1,1].imshow(np.absolute(image - model), vmin=0, vmax=1, interpolation='nearest',origin=0,cmap=cmap,extent=[scalex,-scalex,-scaley, scaley])
    ##### layout and colorbar
    plt.tight_layout()
    #plt.grid()
    plt.subplots_adjust(top=0.85)
    cax = plt.axes([0.93, 0.12, 0.015, 0.7])
    plt.colorbar(im,cax=cax)
    if save!='False':
        fig.savefig(save,dpi=600)
    else:
        plt.show()
    return(None)



##### Main.

if __name__=='__main__':
	print('This is a importable file, see code for more details')
