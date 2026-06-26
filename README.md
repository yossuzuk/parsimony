# Source Code for "Detecting natural selection by parsimony method without inferring ancestral sequences in hemagglutinin of influenza A virus A(H1N1)pdm09"

## Overview
This repository contains C source code used in the study "Detecting natural selection by parsimony method without inferring ancestral sequences in hemagglutinin of influenza A virus A(H1N1)pdm09".

## Directory structure
src/     C source files

## Source files
obsaa.c

sbsnoanc.2.c

sbsnoanc.1.c

## Compilation
gcc adapt-t.log.c -lm -o adapt-t.log
gcc adapt-t.log.negcor.c -lm -o adapt-t.log.negcor
gcc adapt-t.log.oneneg.c -lm -o adapt-t.log.oneneg
gcc adapt-t.log.poscor.c -lm -o adapt-t.log.poscor
gcc ngs.1.c -lm -o ngs.1


## Data
The complete input and output datasets used in this study are available on Zenodo: 

DOI: https://doi.org/

## Contact
Yoshiyuki Suzuki
Nagoya City University

yossuzuk@nsc.nagoya-cu.ac.jp
