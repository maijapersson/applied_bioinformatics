#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# Pipeline for benchmarking
import pybedtools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ruptures as rpt
from remove_artifacts import pred_breakpoints, bed_to_array
from benchmarking import roc_curve, calc_acc_sens_spec, under_estimated_sites, over_estimated_sites

#Extract data
def extract_data(path, row):
    data=pd.read_csv(path, header=None, sep='\t')
    #Extract the individual sample
    individual=np.array(data[row])

    #Extract the reference
    reference=np.array(data[2])
    return individual, reference

#Run the window alogritm
def window(individual, reference, pen_ind, pen_ref, width, min_size):

    #Define the model and algoritm
    model = "l2"  # "l1", "rbf", "linear", "normal", "ar"

    #Fit the data - uncommet this to try different models
    algo = rpt.Window(width=width, model=model, min_size=min_size).fit(individual)
    algo_ref = rpt.Window(width=width, model=model, min_size=min_size).fit(reference)

    #Predict the breakpoints for the individual and reference
    ind_bkps = algo.predict(pen=pen_ind) #- uncommet this to try different penalties etc
    ref_bkps = algo_ref.predict(pen=pen_ref) #- uncommet this to try different models
    return ind_bkps, ref_bkps




# Show the data from ruptures
def display_data(individual, ind_bkps, reference, ref_bkps):
    rpt.display(individual, ind_bkps, figsize=(100, 6))
    plt.title('Before filtering and artifacts removal')
    plt.show()

    #Show the reference
    rpt.display(reference, ref_bkps, figsize=(100, 6))
    plt.title('Reference')
    plt.show()

# calculate and display the filtred result
def display_filtered(ref_bkps , ind_bkps, individual, threshold):
    rem_pred_filter=pred_breakpoints(ref_bkps , ind_bkps, individual, threshold)
    rem_pred_filter_arr=bed_to_array(rem_pred_filter)
    true_bkps_arr=bed_to_array(true_bkps_bed)
    rpt.display(individual, rem_pred_filter_arr, figsize=(100, 6))
    plt.title('After filtering')
    plt.show()
    return rem_pred_filter

#Calculate and display over and under estimated sites.
#Calculate and print the accuracy, sensitivty and specificity
def display_over_under_est(true_bkps_bed, rem_pred_filter, length):
    true_del_with_matches_df = true_bkps_bed.intersect(rem_pred_filter, c=True).to_dataframe()
    false_pos_df=rem_pred_filter.subtract(true_bkps_bed).to_dataframe()
    true_pos_df=rem_pred_filter.intersect(true_bkps_bed).to_dataframe()

    #Under estimated
    under_estimated_sites(true_del_with_matches_df, true_pos_df)
    plt.show()
    #Over estimated
    over_estimated_sites(true_del_with_matches_df, false_pos_df)
    plt.show()
    print("accuracy, sensitivity, specificity")
    print(calc_acc_sens_spec(rem_pred_filter, true_bkps_bed, length))



#1 million 100x
#Define the dataset and individual
individual, reference=extract_data('/Users/linnbeckman/OneDrive/Uppsala Universitet_OneDrive/project_bioinf/Simulated_data/1mil_100x', 5)

#Read the true deletions as a BED file
true_bkps_bed = pybedtools.BedTool("/Users/linnbeckman/OneDrive/Uppsala Universitet_OneDrive/project_bioinf/applied_bioinformatics/1mil_C.bed")


#Calculate break points - example for window method
ind_bkps, ref_bkps=window(individual, reference, 10**3, 10**3, 100, 100)

#Display roc curve and filter
roc_curve(100, ind_bkps, true_bkps_bed, individual, ref_bkps, 1000000)
display_data(individual, ind_bkps, reference, ref_bkps)
rem_pred_filter=display_filtered(ref_bkps , ind_bkps, individual, 24)
display_over_under_est(true_bkps_bed, rem_pred_filter, 1000000)
