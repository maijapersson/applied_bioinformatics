#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Scrip to remove artifacts.
import pandas as pd
import numpy as np
import pybedtools
import matplotlib.pyplot as plt

#Creates a bed file from an array
def create_bed_from_array(array):

    data_array = []
    if len(array) % 2 != 0:
        array=array[:-1] #Remove last value
    for i in range(0, len(array), 2):
        data_array.append(["PA_chr01", array[i],array[i+1]])
    np_array = np.array(data_array)
    dataframe = pd.DataFrame(np_array)
    print()
    return pybedtools.BedTool.from_dataframe(dataframe)

#Removes exact matches from an array
def remove_artifacts_exact(reference, individual):
    ref_set = set(ref)
    deletions = [x for x in ind1 if x not in ref_set]
    return deletions

#Removes matches in a window, w=1000 is the extension of the window.
def remove_artifacts_window(ref_bed, ind_bed):
    ind5_bkps_bed = pybedtools.BedTool(ind_bed)
    ref_bkps_bed = pybedtools.BedTool(ref_bed)
    return ind5_bkps_bed.window(ref_bkps_bed, w=100, v=True)

def bed_to_array(bed):
    # df=pd.read_csv('removed_artifacts_pred.bed',header=None, sep='\t')
    # col_one_arr = df[1].to_numpy()
    # col_two_arr = df[2].to_numpy()
    # start_stop=np.array([0,size])
    # arr = np.concatenate((col_one_arr, col_two_arr))
    # arr = np.concatenate((arr, start_stop))

    # arr=np.sort(arr)

    return pybedtools.BedTool.to_dataframe(bed).loc[:, ["start", "end"]].to_numpy().flatten()

#Filter and remove artifacts
def pred_breakpoints(ref_arr, pred_arr, coverage_array, thr):
        filter_arr=filter_coverage(pred_arr, thr, coverage_array) #Filter some threshold
        pred_filtred_bed=create_bed_from_array(filter_arr) #create bed file
        ref_bed=create_bed_from_array(ref_arr)
        rem_pred_filter=remove_artifacts_window(ref_bed,pred_filtred_bed)
        return rem_pred_filter


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
