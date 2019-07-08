from astropy.io import fits
import gzip
import matplotlib.pyplot as plt
import numpy as np
import re
import math
from scipy.signal import convolve2d
from scipy.interpolate import interp1d
from operator import methodcaller
from FunctionsModule import *

##### Define usefull functions
def openimage(location):
    """
    opens the image and the header
    """
    hdulist = fits.open(location)
    header = hdulist[0].header
    image = np.squeeze(hdulist[0].data)
    return(header, image)

def opendata(filename,parameters,folder):
    """Opens the data from the files
    """
    ##### Unzip and open the data
    # fits
    f_unzip=gzip.open(folder+filename)
    f_unfits=fits.open(f_unzip) 
    datalist=[]
    for i in f_unfits:
        datalist.append(i.data)
    f_unzip.close()
    # params
    par = open(folder+parameters, 'r')
    param = par.readlines()
    par.close()
    return(datalist,param)

def WaveLengthExtract(param,folder):
    """Extract the wavelengths for a .para file
    The function is ugly, but it can be used Blackbox--Like
    """
    # Finding if using custom file
    no_use_other_file = param[9].rstrip('\r\n').split(" ")[4]
    # no custom file
    if no_use_other_file == 'T':
        wavelength_line = param[8].rstrip('\r\n').split(" ")
        n_lambda = float(wavelength_line[2])
        lambda_min = float(wavelength_line[4])
        lambda_max = float(wavelength_line[5])
        # simple linspace
        wavelength = np.exp(np.linspace(np.log(lambda_min),np.log(lambda_max), int(n_lambda))) # because f*** python2
    # custom file
    elif no_use_other_file == 'F':
        file_used_line =  param[10].rstrip('\r\n').split("\t")[0].split(' ')
        file_used = file_used_line[2]
        wavelength_file = open(folder + file_used, 'r')
        wavelength_text = wavelength_file.readlines()
        wavelength_file.close()
        wavelength=[]
        for i in wavelength_text:
            wavelength.append(eval(i[:-1]))
    return(np.array(wavelength))

def plotSED(wl,sed_profile,error=None,label=None,color='g'):
    try :
        error.shape
        plt.errorbar(wl,sed_profile,yerr=2*error,barsabove=True,label=label,fmt='x',ecolor='k',capthick=2.,color=color,alpha=.7,elinewidth=.1)
    except :
        plt.plot(wl,sed_profile,label=label,marker='.',color=color,linewidth=0.5,alpha=.5)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\lambda$ ($\mu$m)')
    plt.ylabel('SED')
    plt.grid(True, which="major", alpha=0.5)
    plt.grid(True, which="minor", alpha=0.1)

def linetodata(line):
    return(re.split(r'\s{2,}',line))

def opendatfile(location):
    f=open(location)
    lines=f.readlines()
    f.close()
    lines=map(linetodata,lines)
    return(list(lines))

def DatafromVOSA(datfile):
    extracted=opendatfile(datfile)
    names=np.array(extracted[7],dtype=str)[1:-1]
    units=np.array(extracted[8],dtype=str)[1:-1]
    data =np.array(extracted[10:],dtype=str)[:,1:-1]
    return(names,units,data)

def SEDfromVOSA(datfile):
    _,_,SEDData=DatafromVOSA(datfile)
    wl,sed_profile,error=SEDData[:,1].astype(float)*1e-4,SEDData[:,4].astype(float),SEDData[:,5].astype(float)
    return(wl,sed_profile,error)

def ergscm2aTOwm2(x,y):
    return(1e1*x*y)

def Gaussian2D(x, y, x_stddev=1.,y_stddev=1.,theta=0., amplitude=1.):
    """
    Redifining the gaussian2D Model of astropy module, for it to be faster in terms of computation. Please refer to the documentation of astropy.modeling.models.Gaussian2D
    """
    c=math.cos(theta)
    s=math.sin(theta)
    c2=math.cos(2.*theta)
    s2=math.sin(2.*theta)
    csq=c**2.
    ssq=s**2.
    a=.5*( csq * (x_stddev )**-2. + ssq * ( y_stddev )**-2.)
    b=.5*s2*(y_stddev**-2.-x_stddev**-2.)
    c=.5*(ssq*(x_stddev**-2.) + csq*(y_stddev**-2.))
    return(amplitude*np.exp( - (a*x**2.+b*(x*y)+c*y**2.) ))

def ExtractBeam(location):
    header,image=openimage(location)
    return(header['BMIN'],header['BMAJ'],header['BPA'])

def ExtractPixtosec(location):
    header,image=openimage(location)
    l1,l2=image.shape
    return(header['CDELT2']*3600,l1)

def MakeBeamMesh(loc1,loc2):
    bmin,bmaj,bpa=ExtractBeam(loc1)
    bmin,bmaj,bpa=bmin*3600/2,bmaj*3600/2,bpa*np.pi/180
    pixtosec,l=ExtractPixtosec(loc2)
    x=np.arange(l)-l/2
    y=np.arange(l)-l/2
    x*=pixtosec
    y*=pixtosec
    x,y=np.meshgrid(x,y)
    beam=Gaussian2D(x,y,bmin,bmaj,bpa,1)
    return(beam)

def Convolve(Im1,Im2):
    return(convolve2d(Im1,Im2))

def CreateFitsFile(location,header,data):
    hdu = fits.PrimaryHDU(data)
    fits.writeto(location,data,header=header,overwrite=True)

def ComputeLimitSignal(image):
    i=image[:500,:500]
    return(3*np.std(i.flatten()))

def ComputeSEDImage(location):
    header,image=openimage(location)
    lim=ComputeLimitSignal(image)
    FluxTot=np.sum(((image>lim)*image).flatten())
    beam_area = np.pi/(4*np.log(2))*header['BMAJ']/header['CDELT2']*header['BMIN']/header['CDELT2']
    Flux=FluxTot/beam_area
    sed=Flux*header['RESTFRQ']*1e-26
    err=ComputeLimitSignal(image)/beam_area*h['RESTFRQ']*1e-26*np.sum(image>lim)/3
    return(sed,err)

def ExperimentalSED(lsed):
    wl=[]
    err=[]
    sed=[]
    origin=[]
    for i in lsed:
        j=np.load(i,allow_pickle=True).item()
        wl.append(j['wl'])
        err.append(j['err'])
        origin.append(j['origin'])
        sed.append(j['sed'])
    return(wl,sed,err,origin)

def Av_extinction(Av, wavelength,extincfile = 'extinction_law.dat'):
    """applies extinction from a file.
    """
    ##### Open Data
    file = open(extincfile, 'r')
    lignes = file.readlines()
    file.close()
    data=np.array([np.array(x.split()[:6],dtype=float) for x in lignes[19:]])
    
    ##### Get absorption
    wl3=data[:,0]
    absorption=data[:,5]
    
    ##### Normalize absoption by its value at 0.547 micron
    def find_closest(val,arr):
        dist=np.abs(arr-val)
        return(np.argmin(dist))
    norm = absorption[find_closest(5.47e-1,wl3)]
    absorption/=norm
    
    ##### Compute the interpolated absorption
    interpoled_absorption=interp1d(wl3,absorption)
    
    ##### compute the new Av
    tau_v = 0.4 * np.log(10) * Av
    correct_Av = np.exp(-tau_v * interpoled_absorption(wavelength))
    return(correct_Av)

def extinct(Av, wavelength, SED):
    return(Av_extinction(Av, wavelength)*SED)


def NormalizedRadialProfile(location,inc,pa,binwidth,par=None,CenterAndWidth=None,save=None):
    inc,pa=float(inc)%360.,float(pa)%360.
    header,image=openimage(location)
    if not par:
        try:
            par=header["PAR"]
            parbool=True
        except:
            par=1000
            parbool=False
    pixtosec=np.abs(3600*header['CDELT1'])
    arcsectoau=1000/par
    pixtoau=arcsectoau*pixtosec
    if CenterAndWidth :
        lx,ly,w=eval(CenterAndWidth)
        image=resizedimage(image,lx,ly,w//2)
    llx,lly=image.shape
    x = np.arange(llx)-llx/2
    y = np.arange(lly)-llx/2
    x, y = np.meshgrid(x, y)
    pixtosec=np.abs(3600*header['CDELT1'])
    xx,yy,rr,angles=deprojectedcoordinates(x,y,inc*degtorad,pa*degtorad)
    m,s=radialbin(image,rr,binwidth)
    n=int(np.max(rr)//binwidth)
    r=np.array(list(range(n+1)))*binwidth
    r*=pixtoau
    m_max=np.nanmax(m)
    m/=m_max
    s/=m_max
    return(r,m,s)

def PlotRadialProfile(r,m,s,unit='normalized',color='k',label=None):
    plt.errorbar(r,m,yerr=s,barsabove=True,linestyle='',markeredgewidth=2,elinewidth=0.5,c=color)
    plt.scatter(r,m,s=0.8,c=color,label=label)
    plt.grid(True)
    plt.title("Average intensity on fixed radius",size=12)
    plt.xlabel('Radius (au)')
    plt.ylabel('Mean flux ({})'.format(unit))
    #plt.yscale('log')
    plt.tight_layout()
    #plt.show()
    #plt.savefig('radialprofile{}.png'.format(name if (name) else int(time())),dpi=600)
    return(None)
