# FCPC_MLs_EFA_SEM
[![DOI](https://zenodo.org/badge/892327534.svg)](https://doi.org/10.5281/zenodo.14207434)
#The following commands and their data is contained :
1. MLs (ELA,AA,RF, and XGBoost)
2. CA (correlation network)
3. EFA
4. SEM (including Shapiro-Wilk test)
5. BayesLiNGAM

#Individual Explanation
1. MLs (ELA,AA,RF, and XGBoost)

   Commands (R) : MLs.R (classification for 8 groups) / MLs_9groups_for_RF_and_XGBoost.R (classification for 9 groups containing Pig_ThB)
   
   Commands (python) : Bubblechart.py　(classification for 8 groups) / Bubblechart_9groups.py (classification for 9 groups containing Pig_ThB) 
   
   Raw data for ELA: ELA_raw.csv; ELA_grouplist.csv
   
   Raw data for AA: AA_binarized_raw_data.csv
   
   Raw data for RF and XGBoost: RF_XG_raw_data.csv (classification for 8 groups) / RF_XG_raw_data_9groups.csv(classification for 9 groups containing Pig_ThB)
   
   Raw data for Bubblechart.py: ML_mix.xlsx (Zip-stored) (original data only)(classification for 8 groups) / ML_mix_new.xlsx (Zip-stored) (new file) (sheet_name="ML_mix_9groups", classification for 9 groups)
   
3. CA (correlation network)

   Commands: CA.R
   
   Raw data for CA: FCPC_MLs.txt

5. EFA

   Commands: EFA.R
   
   Raw data for EFA: FCPC_MLs_AA.csv

7. SEM (including Shapiro-Wilk test)

   Commands: SEM.R

   Raw data for SEM: FCPC_MLs_AA.csv

9. BayesLiNGAM

   Commands: BayesLiNGAM.R

   Raw data for SEM: FCPC_Bayes.csv


 # Please check the license file on this website.
