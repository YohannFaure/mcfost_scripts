from mcfostModule import *
import subprocess

from scipy.interpolate import interp1d
import scipy.optimize as opt

script= """
para=optipar.para
lambda=867

mcfost=~/mcfost_utils/mcfost
rm -rf data_*_old
$mcfost $para
python plotSED.py data_th $para --XPSED --Av 0.6
"""

chmodex='chmod +x {}'

def runmcfost(param):
    filename='optimization.sh'
    f=open(filename,"w")
    f.write(script.format(param,param))
    f.close()
    subprocess.call(chmodex.format(filename),shell=True)
    subprocess.call('./{}'.format(filename),shell=True)
    return(None)

def chi2(inter,x2,y2):
    return(np.sum((1-y2/inter(x2))**2))

lsed=['SEDPoints/SED_J1615_dict.npy','SEDPoints/SED_KOOISTRA_dict.npy','SEDPoints/SED_MAREL_dict.npy']
captions=[r'Our Alma observation',r'Kooistra et al. 2016',r'Van der Marel et al. 2016']

wl,sed,err,origin=ExperimentalSED(lsed)

def cleanersed(listarray):
    t=[]
    for i in listarray:
        if i.shape==():
            t.append(i)
        else:
            for j in i:
                t.append(j)
    return(np.array(t))

x=cleanersed(wl)
y=cleanersed(sed)

def getsed(params):
    ##### Define the working folder
    folder='data_th/'
    filename='sed_rt.fits.gz'
    parameters=params
     ##### Extract the data
    datalist,param=opendata(filename,parameters,folder)
    sed_profile = np.squeeze(datalist[0].data)
    wl=WaveLengthExtract(param,folder)
    sed_profile = extinct(0.6,wl,sed_profile)
    inter = interp1d(wl,sed_profile)
    return(inter)

f=open("empty.para")
emptypar=f.read()
f.close()

def parameterformat(p):
    plist=list(p)
    temp=plist+2*plist[-2:]
    return(emptypar.format(*temp))

def writeparam(p):
    f=open('optipar.para','w')
    f.write(parameterformat(p))
    f.close()
    return(None)

p0=np.array([3e-5,1.5,.4,6.,1.15,
             8e-4,1.5,17.,65.,1.15,
             8e-6,10,95,120,1.15,
             8e-6,20,125,160,1.15,
             0.01,10.,0.01,1000])

p_lim=np.array([
    [1e-8,1e-2],
    [0.5,3.],
    [0.2,2.],
    [2.,10.],
    [1.,1.3],
    [1e-7,1e-2],
    [1.,10.],
    [15.,30.],
    [40.,80.],
    [1.,1.3],
    [1e-7,1e-3],
    [1.,15.],
    [80,100],
    [110,121],
    [1.,1.3],
    [1e-7,1e-3],
    [1.,25.],
    [122,135],
    [136,180],
    [1.,1.3],
    [0.005,1.],
    [1.,1000.],
    [0.005,10],
    [10.,3000]
    ])    

def f_to_minimize(plog):
    p=np.exp(plog)
    writeparam(p)
    runmcfost('optipar.para')
    inter=getsed('optipar.para')
    chi=(chi2(inter,x,y))
    print('\n\n\n',p,'\n\n\n',chi,'\n\n\n')
    return(np.log(chi))

#xxx=opt.minimize(f_to_minimize,np.log(p0),method='TNC',bounds=np.log(p_lim))

def cleanandrandomize(p0list,resultlist,randval=1e-1):
    p0list=list(p0list)
    resultlist=list(resultlist)
    p0listbis=[]
    n=0
    while n<5 and len(p0list)>0:
        amin=np.argmin(resultlist)
        _=resultlist.pop(amin)
        p=p0list.pop(amin)
        if True:#(p<p_lim[:,1]).all() and (p>p_lim[:,0]).all():
            p0listbis.append(p)
            n+=1
    #if n==0:
    #    
    while len(p0listbis)<100:
        p0listbis.append(p0listbis[len(p0listbis)%n])
    return(np.array(p0listbis)*(1+(.5-np.random.random(len(p0list[0])))*randval))

p0list=([p0*(1+(.5-np.random.random(24))*2e-1) for _ in range(100)])

#p0list=np.load('p0list.npy',allow_pickle=True)
#resultlist=np.load('resultlist.npy',allow_pickle=True)

p0list=cleanandrandomize(p0list,resultlist,0.3)

for i in range(5):
    resultlist=[f_to_minimize(np.log(p0list[i])) for i in range(len(p0list))]
    np.save('result_walker_{}.npy'.format(i+5),p0list,allow_pickle=True)
    np.save('result_values_{}.npy'.format(i+5),resultlist,allow_pickle=True)
    p0list=cleanandrandomize(p0list,resultlist)


resultlist=[f_to_minimize(np.log(p0list[i])) for i in range(len(p0list))]
np.save('resultlist2.npy',resultlist,allow_pickle=True)
np.save('p0list2.npy',p0list,allow_pickle=True)
