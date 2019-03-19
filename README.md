# Ribozyme-SeqFun

This repository contain scripts to analyze sequence and structural motifs.

compileRS_FS.m takes all RNAseq and FACSseq analysis data and converts sequences to one-hot vectors, and standardizes activity measures. 

analyzemotifs.m takes the compiled and processed data from compileRS_FS.m and performs pairwise contribution and mutual information analyses using the MatLab scripts findEnt.m, findMutualInformation.m, and findPairwiseMu.m.

ribozyme_XGBR.ipynb is an iPython notebook with the script to perform a XGBoost regression model and compute SHAP values of impact of each feature on the model. 

gbm_RS.R is an R script that runs H2O's AutoML regression model to denoise RNA-Seq and FACS-Seq data. 

genSeqLogos.R contain the script used to generate sequence logos.
