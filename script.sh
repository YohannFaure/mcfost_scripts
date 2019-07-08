#!/bin/bash
# This script can be used to :
# 1) compute the MCFOST simulation of a .para file, specified in line 8
# 2) compute the MCFOST image at lambda fixed in line 9
# 3) plot the SED
# 4) plot the image using FitsSlider (https://github.com/YohannFaure/FitsSlider.git)
# 5) plot the radial profil compared to an original image
para=test.para
lambda=867

### If you are using python instead of python3, change it.
# Set the location of MCFOST
mcfost=~/mcfost_utils/mcfost
# clean folder
rm -rf data_*_old
# Compute mcfost
$mcfost $para
# Plot SED
python3 plotSED.py data_th $para --XPSED --Av 0.6

# compute image
$mcfost $para -img $lambda
# unzip image
gunzip data_$lambda/*.gz
# beam-convolution
python3 ConvolveBeam.py ../DiskFitting/J1615_edit.fits data_$lambda/RT.fits result.fits
# radial plot
python3 DualRadialPlot.py result.fits ../DiskFitting/J1615_edit.fits 45.9835388 326.2477693  6.341708607533389 3. --CenterAndWidth1 '(1555,1555,600)' --CenterAndWidth2 '(1500,1500,600)' --label1 'MCFOST Simulation' --label2 'Alma Data' &
# image show
FitsSlider result.fits
