# Source Code for "Detecting natural selection by parsimony method without inferring ancestral sequences in hemagglutinin of influenza A virus A(H1N1)pdm09"

## Overview

This repository contains C source code used in the study "Detecting natural selection by parsimony method without inferring ancestral sequences in hemagglutinin of influenza A virus A(H1N1)pdm09".

## Directory structure

src/     C source files

## Source files

adapt-t.log.c

adapt-t.log.oneneg.c

ngssitefreq.1.c

ngssitefreq162.1.c

noredunfafsta.1.c

nwktooyako.2.c

obsaa.thresh.1.c

sbsnoanc.1-1.allb.c

sbsnoanc.negcon.1-1.allb.c

sbsnoanc.poscon.1-1.allb.c

takeseq.fas.1.c

transframe1.1.c

## Compilation

gcc adapt-t.log.c -lm -o adapt-t.log

gcc adapt-t.log.oneneg.c -lm -o adapt-t.log.oneneg

gcc ngssitefreq.1.c -lm -o ngssitefreq.1

gcc ngssitefreq162.1.c -lm -o ngssitefreq162.1

gcc noredunfafsta.1.c -lm -o noredunfasta.1

gcc nwktooyako.2.c -lm -o nwktooyako.2

gcc obsaa.thresh.1.c -lm -o obsaa.thresh.1

gcc sbsnoanc.1-1.allb.c -lm sbsnoanc.2-1.c -lm -o sbsnoanc.1.allb

gcc sbsnoanc.negcon.1-1.allb.c -lm sbsnoanc.2-1.c -lm -o sbsnoanc.negcon.1

gcc sbsnoanc.poscon.1-1.allb.c -lm sbsnoanc.2-1.c -lm -o sbsnoanc.poscon.1

gcc takeseq.fas.1.c -lm -o takeseq.fas.1

gcc transframe1.1.c -lm -o transframe1.1

## Data

The complete input and output datasets used in this study are available on Zenodo: 

DOI: https://doi.org/

## Contact

Yoshiyuki Suzuki

Nagoya City University

yossuzuk@nsc.nagoya-cu.ac.jp
