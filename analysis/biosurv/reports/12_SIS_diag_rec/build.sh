#!/bin/sh
qsub -pe smp 16 -cwd -N b12_SIS_diag_rec -j y -b y -l h_vmem=20G make
