#!/bin/sh
#qsub -pe smp 1 -cwd -N r16_SIS_NMF_all -j y -b y -l h_vmem=400G,mem_requested=400G make
qsub -pe smp 1 -cwd -N r16_SIS_NMF_all -j y -b y -l h_vmem=100G,mem_requested=100G make
