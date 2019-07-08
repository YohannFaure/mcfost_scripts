#!/bin/bash

para=test.para
lambda=867

mcfost=~/mcfost_utils/mcfost
rm -rf data_*_old
$mcfost $para
python3 plotSED.py data_th $para --XPSED --Av 0.6

$mcfost $para -img $lambda
gunzip data_$lambda/*.gz
python3 ConvolveBeam.py ../DiskFitting/J1615_edit.fits data_$lambda/RT.fits result.fits
python3 DualRadialPlot.py result.fits ../DiskFitting/J1615_edit.fits 45.9835388 326.2477693  6.341708607533389 3. --CenterAndWidth1 '(1555,1555,600)' --CenterAndWidth2 '(1500,1500,600)' --label1 'MCFOST Simulation' --label2 'Alma Data' &
FitsSlider result.fits

#python3 substractfits.py ../DiskFitting/J1615_edit.fits result.fits resultdiff.fits --norm
