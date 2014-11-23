#!/bin/sh
qsub -pe smp 8 -cwd -N b12_SIS_diag_rec -j y -b y -l h_vmem=60G,mem_requested=60G make -f ../common/SIS_NMF_Makefile DATA_SOURCE=diag_rec
