#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#Scrip to remove artifacts.
import pandas as pd
import pybedtools

#Run the algoritm on reference, the changepoints there should be artifacts. Compare with the individuals and remove matches.
ind1=[435,
 19185,
 19920,
 23580,
 27535,
 96965,
 97405,
 123335,
 129635,
 143400,
 146625,
 243305,
 244185,
 292525,
 293110,
 300000]

ref=[435, 19185, 27535, 143400, 145310, 292380, 293110, 300000]
#Creates a bed file from an array
def create_bed_from_array(arr, name):
    h,t=1,-1;
    arr_noHead_noTail=arr[slice(h,t)] #Remove first and last value
    df_arr=pd.DataFrame(arr_noHead_noTail) #Convert to df
    start=df_arr.iloc[::2, :] #Extract start values
    stop=df_arr.iloc[1::2, :] #Extract stop values
    start.insert(1, "stop", stop) #Merge start & stop
    start=start.rename(index=lambda s: "PA_chr01_1") #Hardcoded row name for testing...
    start.to_csv(name+".bed", header=False,sep='\t') #create bed file


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


create_bed_from_array(ind1, "ind_bkps")
create_bed_from_array(ref, "ref_bkps")


def pred_breakpoints(ref_arr, ind_arr):
    create_bed_from_array(ind_arr, "ind_bkps")
    create_bed_from_array(ref_arr, "ref_bkps")
    remove_artifacts_window( "ref_bkps.bed","ind_bkps.bed").saveas("removed_artifacts_pred.bed")


#remove_artifacts_window( "ref_bkps.bed","ind_bkps.bed").saveas("removed_artifacts_pred.bed")
# print(remove_artifacts_exact(ref,ind1))

pred_breakpoints(ref, ind1)
