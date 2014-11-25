#!/bin/sh
qsub -pe smp 8 -cwd -N r16_SIS_NMF_all -j y -b y -l h_vmem=60G,mem_requested=60G make
