#!/bin/bash
    mcfost=~/mcfost_utils/mcfost
    rm -rf data_*_old
    $mcfost optipar.para
    python3 plotSED.py data_th optipar.para --XPSED
    