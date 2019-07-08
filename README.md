# Some MCFOST Python utilities

This repository contains many python function and bash scripts to use MCFOST properly.

## Installation

### MCFOST
The [Mac](http://ipag.osug.fr/public/pintec/mcfost/mac/mcfost.tgz) and [Linux](http://ipag.osug.fr/public/pintec/mcfost/linux/mcfost.tgz) binaries are here.

Download the `mcfost` file, then :

```
cd path/to/mcfost
export MCFOST_UTILS=~/mcfost_utils/
./mcfost -setup
mv mcfost $MCFOST_UTILS
```

The first time you run MCFOST, make sure to accept the User License and conditions.

### FitSlider

To plot your fits images, I recommand using [FitSlider](https://github.com/YohannFaure/FitsSlider.git).

## The utilities

### ConvolveBeam

Allow you to convolve a image with a beam. It is used to convolve the MCFOST image with the initial image beam.

`python3 ConvolveBeam.py experimental/image/location simulation/image/location output/image/location`

`python3 ConvolveBeam.py J1615_edit.fits data_867/RT.fits result.fits
`

### DualRadialPlot

If you want to plot simultaneously two radial plots, to compare the simulated image and the real image, you can use this script.

It requires many arguments, and you should read it to understand it. In short, it takes the image, the inc,pa,par, the binning width for the profile, and other parameters if you want.

### empty.para

This is super important : it is a para file containing voids, and therefore formatable by python.

To change a para file from python, replace the parameters you want to change by the `{}` string, just load it and read it as any file in Python, and format it. Some functions are built to do it in `mcfostModule.py` and `optimization_mcfost.py`

A rapid string formating tutorial :
```python
# Define a string to format:
>>> teststring='two plus two equal {}, minus one equal {}'

# Formating it, no matter the type of the formater
>>> teststring.format('four','three')
'two plus two equal four, minus one equal three'

>>> teststring.format(4,3)
'two plus two equal 4, minus one equal 3'

# Formating it with an iterable (tuple, list, etc.) takes one more step : a star
>>> iterable=(4,3)
>>> teststring.format(iterable)
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
IndexError: tuple index out of range

>>> teststring.format(*iterable)
'two plus two equal 4, minus one equal 3'
```

### extinction_law.dat

Used to compute the extinction.

### FunctionsModule.py

Comes from [DiskFitting](https://github.com/YohannFaure/DiskFitting.git).

### launch.sh

Just a script to launch the `optimization_mcfost.py` on Leftraru.

### mcfostModule.py

Contains most of the useful functions of that project.

### optimization_mcfost.py

An example of usage of the python functions, used to optimize MCFOST parameters to fit an experimental SED.

### optimization.sh and optipar.para

Output files of `optimization_mcfost.py`.

### Paramfiles.py

Complex (and legacy) parameters manipulation program, found online.

### ParamFinal.para

The best parameters I found during my internship. It still needs some workaround, especially in the µm zone of the SED.

### PlotSED.py

Just an SED plotting utility.

`python3 plotSED.py data_th params.para --XPSED --Av 0.6`

### README.md

This file

### script.sh

This script can be used to :

* compute the MCFOST simulation of a .para file, specified in line 8
* compute the MCFOST image at lambda fixed in line 9
* plot the SED
* plot the image using [FitsSlider](https://github.com/YohannFaure/FitsSlider.git)
* plot the radial profile compared to an original image

### test.para

A parameters file.

### vosaresults.dat

The SED as given by the VOSA database.

### SEDPoints

Contains dictionnary, saved using numpy, with the different SED points found in the litterature

Format of the dictionnaries :
```python
>>> x=np.load('SED_J1615_dict.npy')
>>> dic=x.item()
>>> dic.keys()
dict_keys(['wl', 'err', 'sed', 'origin'])
```
* wl: the wavelengths
* sed: the SED points
* err: the errors on the SED points
* origin: a string describing the data (example: doi reference)
