#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 12:48:31 2021

@author: annabelbeichman
"""
from SigProfilerExtractor import sigpro as sig
from datetime import date
import os

todaysdate = date.today().strftime('%Y%m%d')

print("Today's date:", todaysdate)
# matrix file in sigprofiler format (with central A mutation types rev comped to T)
reference_genome='GRCh37' # either mouse or human are the only ones supported. run with both? supposed to not matter but still says the genome in the
# run information -- GRCh37 or mm10 

matrixFile="/Users/annabelbeichman/Documents/UW/BearProject/results/mutyper/20200916/mutyperResults_20210401_3mer_nostrict/mapped_to_brown_bear/summaryTables/All_Intervals_kmerCountsForSigProfilerExtractor.ALLFREQUENCIESSEPARATED.bears.txt"
# https://osf.io/t6j7u/wiki/2.%20Quick%20Start%20Example/
label="bears_kSFSAllFreqsSep_"+reference_genome+"GenomeRef"
outdir="/Users/annabelbeichman/Documents/UW/DNAShapeProject/results/SigProfilerExtractor/"+str(todaysdate)+"_SigProfilerExtractorResults_"+label+"/"


if not os.path.exists(outdir):
    os.makedirs(outdir) # for saving plots; need mkdirs for it to be recursive



# so the docs say that ref genome only matters if working from vcf
# but there are also reference signatures that differ for each genome (mm10 vs human)
# so how do those come into play?
# trying to run both ways and see what happens? are there different reference signatures?
# if you leave it blank then human is default
# I am getting different answers with human vs mouse ! tell Kelley for her work 
# upping replicates to 500 from 100 to hopefully get more stability
sig.sigProfilerExtractor(input_type="matrix",output=outdir, input_data=matrixFile, reference_genome= reference_genome, minimum_signatures=1, maximum_signatures=10, nmf_replicates=500, cpu=-1)
# trying with mm10 instead 
# seems like ref genome doesn't matter but run with both genomes just in case