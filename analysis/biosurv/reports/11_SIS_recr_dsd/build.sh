#!/bin/sh
qsub -pe smp 16 -cwd -N b11_SIS_recr_dsd -j y -b y -l h_vmem=20G make
