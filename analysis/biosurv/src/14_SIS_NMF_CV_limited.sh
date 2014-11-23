#!/bin/bash
#ulimit -v 268435456
#ulimit -v 260000000
ulimit -v 200000000
/home/marpin/bin/Rscript 14_SIS_NMF_CV.R | tee > 14_SIS_NMF_CV.Rout
