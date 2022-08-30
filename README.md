# MooN-matching
Date: 30.08.2022

This repository contains the R code to implement the direct M-out-of-N bootstrap for nearest neighbor matching along with the code used in the simulation study in:

Walsh, C.P. and Jentsch, C. (2022+), "Nearest neighbor matching: Does direct M-out-of-N bootstrapping without bias correction work when the naive bootstrap fails"


This repository contains:

* R code in the file **moon_matching.r** to implement our M-out-of-N bootstrap. Further details are given in the file.
* R code in the file **simulation_function.r** to do one simulation run for a fixed ratio of treated to control units over various sample sizes and bootstrap resample sizes. Further details are given in the file.   
* R code used in the simulation study of the paper for the ratio of treated to control given by *alpha = 0.2, 1, 2.84* and *5* is in the files **main_file_alpha_02.r**, **main_file_alpha_1.r**, **main_file_alpha_bar.r**, **main_file_alpha_5.r** respectively. 
* the estimation results from the simulation study are contained in the folder **Results/** in the files **ATET_MooN_alpha_0.2.RData**, **ATET_MooN_alpha_1.RData**, **ATET_MooN_alpha_2.84.RData** and **ATET_MooN_alpha_5.RData**  
* the folder **Results/** also contains the file **summarize_results_paper.r** to construct:
   * **paper_plot_point_estimates.pdf**, which corresponds to Figrue 1 of the paper  
   * and **Coverage_Table.tex**, which gives the coverage probabilities in Table 1 of the paper.   


ACKNOWLEDGEMENTS: The simulations were run on the supercomputer CLAIX at RWTH Aachen University as part of the NHR4CES infrastructure with computing resources under the project p0020105.

CAUTION: The simulations in the files **main_file_alpha_02.r**, **main_file_alpha_1.r**, **main_file_alpha_bar.r**, **main_file_alpha_5.r** were run on 48 cores with run times in hours of 07:43:36, 14:00:32,  10:59:47 and 07:56:22 respectively.   






