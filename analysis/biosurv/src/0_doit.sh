#!/bin/sh
test ! -e 01_cpv_prep.done && Rscript 01_cpv_prep.R > 01_cpv_prep.Rout && touch 01_cpv_prep.done
test ! -e 02_gex_process.done && Rscript 02_gex_process.R > 02_gex_process.Rout && touch 02_gex_process.done
test ! -e 03_subset_gex_cpv.done && Rscript 03_subset_gex_cpv.R > 03_subset_gex_cpv.Rout && touch 03_subset_gex_cpv.done
test ! -e 04_norm_gex.done && Rscript 04_norm_gex.R > 04_norm_gex.Rout && touch 04_norm_gex.done
test ! -e 05_match_gex_cpv.done && Rscript 05_match_gex_cpv.R > 05_match_gex_cpv.Rout && touch 05_match_gex_cpv.done
test ! -e 06_EDA.done && Rscript 06_EDA.R > 06_EDA.Rout && touch 06_EDA.done
test ! -e 07_SIS_prep.done && Rscript 07_SIS_prep.R > 07_SIS_prep.Rout && touch 07_SIS_prep.done
