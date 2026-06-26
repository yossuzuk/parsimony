# Source code for "Detecting natural selection by parsimony method without inferring ancestral sequences in hemagglutinin of influenza A virus A(H1N1)pdm09"

## Overview

This repository contains C source code used in the study "Detecting natural selection by parsimony method without inferring ancestral sequences in hemagglutinin of influenza A virus A(H1N1)pdm09".

## Directory structure

src/     C source files

## Source files

adapt-t.log.c : conducting two-tailed test for natural selection

adapt-t.log.negcor.c : conducting correction for multiple testing of negative selection

adapt-t.log.oneneg.c : conducting one-tailed test for negative selection

adapt-t.log.poscor.c : conducting correction for multiple testing of positive selection

ancseqpamp.1.c : inferring ancestral sequence at each interior node using pamp in paml

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

sbsstr.1-5.allb.c : computing ss, sn, cs, and cn at each codon site with inferring ancesltral sequences

sbsstr.2-5.c : computing ss, sn, cs, and cn at each codon site with inferring ancesltral sequences

sbsstr.h : computing ss, sn, cs, and cn at each codon site with inferring ancesltral sequences

takeseq.fas.1.c : changing the format of alignment

transframe1.1.c : translating nucleotide sequences into amino acid sequences based on standard genetic code

## Compilation

gcc adapt-t.log.c -lm -o adapt-t.log

gcc adapt-t.log.negcor.c -lm -o adapt-t.log.negcor

gcc adapt-t.log.oneneg.c -lm -o adapt-t.log.oneneg

gcc adapt-t.log.poscor.c -lm -o adapt-t.log.poscor

gcc ancseqpamp.1.c -lm -o ancseqpamp.1

gcc ngssitefreq.1.c -lm -o ngssitefreq.1

gcc ngssitefreq162.1.c -lm -o ngssitefreq162.1

gcc noredunfafsta.1.c -lm -o noredunfasta.1

gcc nwktooyako.2.c -lm -o nwktooyako.2

gcc obsaa.thresh.1.c -lm -o obsaa.thresh.1

gcc sbsnoanc.1-1.allb.c -lm sbsnoanc.2-1.c -lm -o sbsnoanc.1.allb

gcc sbsnoanc.negcon.1-1.allb.c -lm sbsnoanc.2-1.c -lm -o sbsnoanc.negcon.1.allb

gcc sbsnoanc.poscon.1-1.allb.c -lm sbsnoanc.2-1.c -lm -o sbsnoanc.poscon.1.allb

gcc sbsstr.1-5.allb.c -lm sbsstr.2-5.c -lm -o sbsstr.5.allb

gcc takeseq.fas.1.c -lm -o takeseq.fas.1

gcc transframe1.1.c -lm -o transframe1.1

## Example

adapt-t.log h1.19949.sbsnoanc.1.allb.ml.-0.1.5.0 > h1.19949.sbsnoanc.1.allb.ml.-0.1.5.0.adt

adapt-t.log.negcor h1.19949.sbsnoanc.1.allb.ml.-0.1.5.0 > h1.19949.sbsnoanc.1.allb.ml.-0.1.5.0.negcor

adapt-t.log.oneneg h1.19949.sbsnoanc.1.allb.ml.0.0.5.0 > h1.19949.sbsnoanc.1.allb.ml.0.0.5.0.adt.on

adapt-t.log.poscor h1.19949.sbsnoanc.1.allb.ml.-0.1.5.0 > h1.19949.sbsnoanc.1.allb.ml.-0.1.5.0.poscor

ancseqpamp.1 h1.19949.gtrgi.ml.acc.oyako seqfile AF112741.1 > h1.gtrgi.ml.acc.ancseqpamp

ngssitefreq.1 2008.orf > 2008.ngs

ngssitefreq162.1 2008.orf > 2008.ngs.162

noredunfasta.1 h1.37433+1.orf > h1.19949+1.orf

nwktooyako.2 h1.19949+1.gtrgi.ml.acc.nwk > h1.19949+1.gtrgi.ml.acc.oyako

obsaa.thresh.1 h1.19949.orf 0.0 > h1.19949.obsaa.0.0

sbsnoanc.1.allb h1.19949.orf.tk h1.19949+1.gtrgi.ml.acc.oyako h1.19949.obsaa.-0.1 5.0 > h1.19949.sbsnoanc.1.allb.ml.-0.1.5.0

sbsnoanc.negcon.1.allb h1.19949.orf.tk h1.19949+1.gtrgi.ml.acc.oyako h1.19949.obsaa.-0.1 5.0 > h1.19949.sbsnoanc.negcon.1.allb.ml.-0.1.5.0

sbsnoanc.poscon.1.allb h1.19949.orf.tk h1.19949+1.gtrgi.ml.acc.oyako h1.19949.obsaa.-0.1 5.0 > h1.19949.sbsnoanc.poscon.1.allb.ml.-0.1.5.0

sbsstr.5.allb h1.gtrgi.ml.acc.ancseqpamp h1.19949+1.gtrgi.ml.acc.oyako h1.19949.obsaa.-0.1 5.0 > h1.19949.sbsstr.5.allb.ml.-0.1.5.0

takeseq.fas.1 h1.19949+1.orf > h1.19949+1.orf.tk

transframe1.1 h1.37433.orf > h1.37433.aa

## Data

The input and output datasets used in this study are available on Zenodo: 

DOI: https://doi.org/10.5281/zenodo.20932236


## Contact

Yoshiyuki Suzuki

Nagoya City University

yossuzuk@nsc.nagoya-cu.ac.jp
