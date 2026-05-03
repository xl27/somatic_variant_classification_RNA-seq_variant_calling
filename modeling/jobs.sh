#!/bin/bash
#SBATCH --job-name=weights_Filter-tcga_perDisease
#SBATCH --mem=100G 
#SBATCH --time=1-00:00               
#SBATCH --output=xgboost-weights_Filter-tcga_perDisease.log 


module load Anaconda3
source activate ml-env
conda activate ml-env

#python xgboost-weights_noSampleFilter.py
#python xgboost-SMOTE_noSampleFilter.py
#python xgboost-weights_Filter-num20.py
#python xgboost-SMOTE_Filter-num20.py
#python xgboost-weights_Filter-tcga.py
#python xgboost-SMOTE_Filter-tcga.py
python xgboost-weights_Filter-tcga_perDisease.py
#python xgboost-SMOTE_Filter-tcga_perDisease.py