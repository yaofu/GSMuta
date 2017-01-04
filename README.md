############# libraries for R
library(abind)
library(foreach)
library(doParallel)
library(parallel)
library("data.table")
library(poibin)
library(survcomp)
library(methods)
library(MASS)



############ Code #############
the underlying codes


########### Data ################
data required for the code



########### Test data ###############
Two datasets:
	TCGA_BRCA.maf
	TCGA_GBM.maf 





############ Running Script ##########
GSMuta.sh



Examples: 
bash GSMuta.sh       (To see how to use it)
bash GSMuta.sh -i test_data/TCGA_BRCA.maf      (General run)
bash GSMuta.sh -i test_data/TCGA_GBM.maf -bed /net/kodiak/volumes/delta/shared/home/yao_new/MADGiC_mod/GBM_coverage -core 2  (Test with coverage. Caution: pls set 'core' to 2 (small values), otherwise it will screw up the memory...)  
