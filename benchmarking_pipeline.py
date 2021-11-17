#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# Pipeline for benchmarking
import pybedtools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ruptures as rpt
from remove_artifacts import pred_breakpoints
from benchmarking import roc_curve, calc_acc_sens_spec, under_estimated_sites, over_estimated_sites

#Define the dataset and wich individual
x100_cov_1mil = pd.read_csv('/Users/linnbeckman/OneDrive/Uppsala Universitet_OneDrive/project_bioinf/Simulated_data/1mil_100x', header=None, sep='\t')
individual_row=5

#Extract the individual sample
individual=np.array(x100_cov_1mil[individual_row])

#Extract the reference
reference=np.array(x100_cov_1mil[2])

#Read the true deletions as a BED file
true_bkps_bed = pybedtools.BedTool("1mil_C.bed")

#Define the model and algoritm
model = "l2"  # "l1", "rbf", "linear", "normal", "ar"
min_size=50 #Minimum size between the deletions

#Fit the data
# algo = rpt.BottomUp(model=model, min_size=min_size).fit(individual)
# algo_ref = rpt.BottomUp(model=model, min_size=min_size).fit(reference)

#Predict the breakpoints for the individual and reference
# ind_bkps = algo.predict(pen=2000000)
ind_bkps=[12320, 16415, 45830, 46935, 190425, 191035, 289420, 290640, 363215, 365350, 434745, 442560, 510735, 515010, 520990, 527215, 735955, 745230, 854485, 860650, 884760, 886285, 932730, 935845, 1000000]
# #
# ref_bkps = algo_ref.predict(pen=2000000)
ref_bkps=[510735, 515010, 1000000]

#ROC curve with different filters
roc_curve(100, ind_bkps, true_bkps_bed, individual, ref_bkps, 1000000)


rem_pred_filter=pred_breakpoints(ref_bkps , ind_bkps, individual, 1)
true_del_with_matches_df = true_bkps_bed.intersect(rem_pred_filter, c=True).to_dataframe()
false_pos_df=rem_pred_filter.subtract(true_bkps_bed).to_dataframe()
true_pos_df=rem_pred_filter.intersect(true_bkps_bed).to_dataframe()
#Under estimated
under_estimated_sites(true_del_with_matches_df, true_pos_df)

#Over estimated
over_estimated_sites(true_del_with_matches_df, false_pos_df)

#Accuracy, sensitivity and specificity for choosen parameters
