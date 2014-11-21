#!/bin/sh
qsub -pe smp 16 -cwd -N b09_SIS_diag_dsd -j y -b y make
