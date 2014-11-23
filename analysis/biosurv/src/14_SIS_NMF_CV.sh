#!/bin/sh
qsub -pe smp 8 -N s14_SIS_NMF_CV -cwd -j y -b y -l h_vmem=60G,mem_requested=60G /home/marpin/bin/Rscript 14_SIS_NMF_CV.R
