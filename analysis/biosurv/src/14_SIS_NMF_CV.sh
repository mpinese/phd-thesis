#!/bin/sh
qsub -pe smp 32 -N s14_SIS_NMF_CV -cwd -j y -b y -l h_vmem=8G,mem_requested=8G ./14_SIS_NMF_CV_limited.sh
