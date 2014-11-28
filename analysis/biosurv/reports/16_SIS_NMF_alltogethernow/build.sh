#!/bin/sh
qsub -pe smp 4 -cwd -N r16_SIS_NMF_all -j y -b y -l h_vmem=50G,mem_requested=50G make
