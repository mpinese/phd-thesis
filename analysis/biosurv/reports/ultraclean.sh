#!/bin/sh
cd 09_SIS_diag_dsd/ && make -f ../common/SIS_NMF_Makefile ultraclean && cd ..
cd 10_SIS_surg_dsd/ && make -f ../common/SIS_NMF_Makefile ultraclean && cd ..
cd 11_SIS_recr_dsd/ && make -f ../common/SIS_NMF_Makefile ultraclean && cd ..
cd 12_SIS_diag_rec/ && make -f ../common/SIS_NMF_Makefile ultraclean && cd ..
cd 13_SIS_surg_rec/ && make -f ../common/SIS_NMF_Makefile ultraclean && cd ..
