#!/bin/sh
qsub -pe smp 16 -cwd -N b13_SIS_surg_rec -j y -b y make
