#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Scrip to remove artifacts.
import pandas as pd
import numpy as np
import pybedtools
from benchmarking import calc_acc_sens_spec
import matplotlib.pyplot as plt

#Creates a bed file from an array
def create_bed_from_array(array):




    data_array = []
    for i in range(0, len(array), 2):
        data_array.append(["PA_Chr01", array[i],array[i+1]])
    np_array = np.array(data_array)
    dataframe = pd.DataFrame(np_array)
    print()
    return pybedtools.BedTool.from_dataframe(dataframe)
    # if len(arr) % 2 == 0:
    #     h,t=1,-1;
    #     arr_noHead_noTail=arr[slice(h,t)] #Remove first and last value
    # else:
    #     arr_noHead_noTail=arr[:-1] #Remove last value

    # df_arr=pd.DataFrame(arr) #Convert to df
    #
    # if len(dataframe)>0:
    #
    #     start=df_arr.loc[::2, :] #Extract start values
    #     stop=df_arr.loc[1::2, :] #Extract stop values
    #     start.insert(1, "stop", stop) #Merge start & stop
    #     start.insert(loc=0, column='chrom', value="PA_chr01")
    #     #start=start.rename(index=lambda s: "PA_chr01") #Hardcoded row name for testing...
    #     #start.to_csv(name+".bed", header=False,sep='\t') #create bed file
    #     #bedfile=pybedtools.BedTool(name+".bed")


#Removes exact matches from an array
def remove_artifacts_exact(reference, individual):
    ref_set = set(ref)
    deletions = [x for x in ind1 if x not in ref_set]
    return deletions

#Removes matches in a window, w=1000 is the extension of the window.
def remove_artifacts_window(ref_bed, ind_bed):
    ind5_bkps_bed = pybedtools.BedTool(ind_bed)
    ref_bkps_bed = pybedtools.BedTool(ref_bed)
    return ind5_bkps_bed.window(ref_bkps_bed, w=1000, v=True)

def bed_to_array(bed, size):
    df=pd.read_csv('removed_artifacts_pred.bed',header=None, sep='\t')
    col_one_arr = df[1].to_numpy()
    col_two_arr = df[2].to_numpy()
    start_stop=np.array([0,size])
    arr = np.concatenate((col_one_arr, col_two_arr))
    arr = np.concatenate((arr, start_stop))
    arr=np.sort(arr)
    return arr

def pred_breakpoints(ref_arr, ind_arr):
    ind=create_bed_from_array(ind_arr)
    ref=create_bed_from_array(ref_arr)
    remove_artifacts_window( ref,ind).saveas("removed_artifacts_pred.bed")
    return bed_to_array("removed_artifacts_pred.bed", 1000000)


def filter_coverage(pred_array, threshold, coverage_array):
    avg_cov_arr=[]
    pred_filterd=[]
    for i in range(len(pred_array)-1):
        start=pred_array[i]
        stop=pred_array[i+1]

        temp=coverage_array[start: stop]

        avg_coverage=sum(temp)/len(temp)
        if avg_coverage<threshold:
            pred_filterd.extend([start, stop])

        avg_cov_arr.append(avg_coverage)
    return pred_filterd


def roc_curve(avg_coverage, pred_array, true_del_bed, coverage_array, ref_arr):

    threshold = np.linspace(0, 0.2*avg_coverage,20)
    acc_arr = []
    spec_arr = []
    sens_arr = []
    for thr in threshold:
        filter_arr=filter_coverage(pred_array, thr, coverage_array) #Filter some threshold

        # pred_bkps=pred_breakpoints(ref_arr, filter_arr) #remove artifacts

        pred_filtred_bed=create_bed_from_array(filter_arr) #create bed file

        accuracy, sensitivity, specificity  = calc_acc_sens_spec(pred_filtred_bed, true_del_bed)
        acc_arr.append(accuracy)
        spec_arr.append(specificity)
        sens_arr.append(sensitivity)

    print("Sensitivity: ")
    print(sens_arr)
    print("specificity")
    print(spec_arr)
    print("accuracy")
    print(acc_arr)
    plt.figure(1)
    plt.plot(1-np.array(spec_arr), sens_arr)
    #plt.plot(threshold, acc_arr)
    plt.show()



arr=[12320, 16475, 289420, 290640, 363280, 365350, 434810, 442620, 520990, 527215, 735955, 745230, 854485, 860590, 932730, 935905]
