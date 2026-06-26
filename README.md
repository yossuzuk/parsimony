# Source Code for "Detecting natural selection by parsimony method without inferring ancestral sequences in hemagglutinin of influenza A virus A(H1N1)pdm09"

## Overview

This repository contains C source code used in the study "Detecting natural selection by parsimony method without inferring ancestral sequences in hemagglutinin of influenza A virus A(H1N1)pdm09".

## Directory structure

src/     C source files

## Source files

adapt-t.log.c : conducting two-tailed test for natural selection

adapt-t.log.oneneg.c : conducting one-tailed test for negative selection

ngssitefreq.1.c : counting n-linked glycosylation sites starting at each position

ngssitefreq162.1.c : counting n-linked glycosylation sites at positions 162-164

noredunfafsta.1.c : eliminating redundant sequences

nwktooyako.2.c : changing the format of phylogenetic tree

obsaa.thresh.1.c : making the list of amino acids observed with proprtion higher than threshold value at each site

sbsnoanc.1-1.allb.c : computing ss, sn, cs, and cn at each codon site without inferring ancesltral sequences

sbsnoanc.2-1.c : computing ss, sn, cs, and cn at each codon site without inferring ancesltral sequences

sbsnoanc.h : computing ss, sn, cs, and cn at each codon site without inferring ancesltral sequences

sbsnoanc.negcon.1-1.allb.c : computing ss, sn, cs, and cn at each codon site without inferring ancesltral sequences for conservative detection of negative selection

sbsnoanc.poscon.1-1.allb.c : computing ss, sn, cs, and cn at each codon site without inferring ancesltral sequences for conservative detection of positive selection

takeseq.fas.1.c : changing the format of alignment

transframe1.1.c : translating nucleotide sequences into amino acid sequences based on standard genetic code

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
