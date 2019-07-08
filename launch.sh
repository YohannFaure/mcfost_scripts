#!/bin/bash
#SBATCH --partition=largemem
#SBATCH --ntasks=44
#SBATCH --ntasks-per-node=44
#SBATCH --mem-per-cpu=4365
#SBATCH --mail-user=faure.yohann@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --output=results/output.out
#SBATCH --error=results/output.err
#SBATCH --job-name=mcfost

export MCFOST_UTILS=~/mcfost_utils/
ml Python/3.7.2
ml matplotlib/3.0.2-Python-3.7.2
cd ~/mcfost

function python3 () {
    python
}

python3 optimization_mcfost.py
