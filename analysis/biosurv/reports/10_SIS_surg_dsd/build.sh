#!/bin/sh
qsub -pe smp 16 -cwd -N b10_SIS_surg_dsd -j y -b y -l h_vmem=20G make
