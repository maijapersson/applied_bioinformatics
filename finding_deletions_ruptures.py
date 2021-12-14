#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import pybedtools
import matplotlib.pyplot as plt
import ruptures as rpt
import time
from sys import argv
#Calculate the average coverage
def calculate_avg_thr(path):
    avg_thr=[]
    data=pd.read_csv(path, header=None, sep='\t')
    #Extract the individual sample
    for i in range(len(data.columns)-2):
        avg_thr.append(3*np.mean(np.array(data[i+2])))
    return avg_thr


#Creates a bed file from an array
def create_bed_from_array(array, chr_name):

    data_array = []
    if len(array) % 2 != 0:
        array=array[:-1] #Remove last value
    for i in range(0, len(array), 2):
        data_array.append([chr_name, array[i],array[i+1]])
    np_array = np.array(data_array)
    dataframe = pd.DataFrame(np_array)
    return pybedtools.BedTool.from_dataframe(dataframe)


#Removes matches with a 60% match
def remove_artifacts(ref_bed, ind_bed):
    ind5_bkps_bed = pybedtools.BedTool(ind_bed)
    ref_bkps_bed = pybedtools.BedTool(ref_bed)
    return ind5_bkps_bed.subtract(ref_bkps_bed, f=0.6, r=True, A=True )#

def bed_to_array(bed):
    # df=pd.read_csv('removed_artifacts_pred.bed',header=None, sep='\t')
    # col_one_arr = df[1].to_numpy()
    # col_two_arr = df[2].to_numpy()
    # start_stop=np.array([0,size])
    # arr = np.concatenate((col_one_arr, col_two_arr))
    # arr = np.concatenate((arr, start_stop))

    # arr=np.sort(arr)
    return pybedtools.BedTool.to_dataframe(bed).loc[:, ["start", "end"]].to_numpy().flatten()

#Calculates and filter the average coverage for every interval
def filter_coverage(pred_array, coverage_array):
    avg_cov_arr=[]
    pred_filterd=[]
    for i in range(len(pred_array)-1):
        start=pred_array[i]
        stop=pred_array[i+1]

        temp=coverage_array[start: stop]

        avg_coverage=sum(temp)/len(temp)
        avg_cov_arr.append(avg_coverage)

    (n, bins, patches) = plt.hist(avg_cov_arr, bins=25)
    #print("threshold: ")
    #print(bins[1])
   # plt.ylabel('Number of intervals')
   # plt.xlabel('Average coverage per interval')
    #plt.show()

    for i in range(len(pred_array)-1):
        start=pred_array[i]
        stop=pred_array[i+1]

        if avg_cov_arr[i]<(bins[1]):
            pred_filterd.extend([start, stop])

    return pred_filterd, bins[1]

#Creates an info array about the predicted change points
def create_info_array(bedObject):
        info=[]
        for feature in bedObject:
            info.append(feature.stop-feature.start)
        min_info=np.min(info)
        max_info=np.max(info)
        avg_info=np.mean(info)
        len_info=len(info)
        return min_info, max_info, avg_info, len_info

#Calculates the changepoints
def predict_deletions(path, chr_name, pen_ref, pen_ind, arr_avg_cov, ref_index, sample_index='all'):


    #Read in the data
    data=pd.read_csv(path, header=None, sep='\t')

    #Define the model and algoritm
    model = "l2"  # "l1", "rbf", "linear", "normal", "ar"
    min_size=100

    #Extract the reference
    all_ref=[]

    #Iterate over all references
    for index in ref_index:

        reference=np.array(data[index])
        reference[reference > arr_avg_cov[index-2]] = arr_avg_cov[index-2]
        t0 = time.time()
        algo_ref = rpt.BottomUp(model=model, min_size=min_size).fit(reference)
        ref_bkps = algo_ref.predict(pen=pen_ref)
        t1 = time.time()
        #print("Reference "+str(index)+": ")
        #print(t1-t0)

        (filter_array_ref, threshold)=filter_coverage(ref_bkps, reference) #Filter the reference
        #Make plot and save it
        rpt.display(reference, filter_array_ref, figsize=(100, 6))
        plt.title('Reference ' + str(index))
        plt.xlabel('Position')
        plt.ylabel('Coverage')
        plt.xticks(np.arange(0, len(reference)+1, 10000))
        plt.savefig('results/reference_'+ str(index) + '.jpg')


        all_ref.extend(filter_array_ref)

        #get info
        ref_bed_1=create_bed_from_array(filter_array_ref, chr_name) #create bed file and save it
        (min_info, max_info, avg_info, len_info)=create_info_array(ref_bed_1)

        #Write to file
        info_file = open("results/info_reference_" + str(index) + ".txt", "w")

        info_file.write("This is a file containing information about reference with index " + str(index) + ". \n")
        info_file.write("Smallest artifact: " + str(min_info) + "\n")
        info_file.write("Biggest artifact: " + str(max_info) + "\n")
        info_file.write("Average artifact: " + str(avg_info)+ "\n")
        info_file.write("Number of artifacts: " + str(len_info)+ "\n")
        info_file.write("Artifact threshold: " + str(threshold)+ "\n")
        info_file.write("Average coverage * 3 cut off: " + str(arr_avg_cov[index-2]) + "\n")
        info_file.close()

    ref_bed=create_bed_from_array(all_ref, chr_name) #create bed file and save it
    ref_bed_sort=ref_bed.sort().merge().saveas("results/reference.bed")

    #Extract the individual sample
    if sample_index=='all': #If 'all' is specifed, the predictions will be done for all columns
        sample_index=range(2, len(data.columns))

    for i in sample_index:

        if i not in ref_index:
            individual=np.array(data[i])
            individual[individual > arr_avg_cov[i-2] ] = arr_avg_cov[i-2]

            #Fit the data
            t0 = time.time()
            algo = rpt.BottomUp(model=model, min_size=min_size).fit(individual)

            #Predict the breakpoints for the individual and reference
            ind_bkps = algo.predict(pen=pen_ind)
            t1 = time.time()
            #print("Sample:" + str(i) + ', running time: ')
            #print(t1-t0)

            fig = plt.figure(num = "ruptures_figure", figsize=(100, 10))
            #define grid of nrows x ncols
            gs = fig.add_gridspec(3, 1)
            _, curr_ax = rpt.display(individual, ind_bkps, num="ruptures_figure")

            curr_ax[0].set_position(gs[0].get_position(fig))
            curr_ax[0].set_subplotspec(gs[0])

            plt.title('Sample '  + str(i)+ ' before filtration and removal of artifacts')
            plt.xlabel('Position')
            plt.ylabel('Coverage')
            plt.xticks(np.arange(0, len(individual)+1, 10000))

            #filter the sample
            (filter_array, threshold)=filter_coverage(ind_bkps, individual)

            _, curr_ax = rpt.display(individual, filter_array,  num="ruptures_figure")
            curr_ax[0].set_position(gs[1].get_position(fig))
            curr_ax[0].set_subplotspec(gs[1])

            plt.title('Sample '  + str(i)+ '  After filtration')
            plt.xlabel('Position')
            plt.ylabel('Coverage')
            plt.xticks(np.arange(0, len(individual)+1, 10000))

            pred_filtred_bed=create_bed_from_array(filter_array, chr_name) #create bed file

           # Remove artifacts from reference
            rem_bed=remove_artifacts(ref_bed_sort,pred_filtred_bed).saveas("results/sample_" + str(i) + ".bed") # save it

            rem_pred_filter_arr=bed_to_array(rem_bed) #Remove reference bkps and convert to array.

            _, curr_ax = rpt.display(individual, rem_pred_filter_arr, num="ruptures_figure")
            curr_ax[0].set_position(gs[2].get_position(fig))
            curr_ax[0].set_subplotspec(gs[2])
            plt.title('Sample '  + str(i)+ '  After filtration and removal of artifacts')
            plt.xlabel('Position')
            plt.ylabel('Coverage')
            plt.xticks(np.arange(0, len(individual)+1, 10000))
            fig.savefig('results/sample_'+ str(i) + '.jpg')


            #print("Coverage cutoff: ")
            #print(arr_avg_cov[i-2])

            (min_info, max_info, avg_info, len_info)=create_info_array(rem_bed)

            #Write to file
            info_file = open("results/info_sample_" + str(i) + ".txt", "w")

            info_file.write("This is a file containing information about sample with index " + str(i) + ". \n")
            info_file.write("Smallest deletion: " + str(min_info) + "\n")
            info_file.write("Biggest deletion: " + str(max_info) + "\n")
            info_file.write("Average deletion: " + str(avg_info)+ "\n")
            info_file.write("Number of deletions: " + str(len_info)+ "\n")
            info_file.write("Deletion threshold: " + str(threshold)+ "\n")
            info_file.write("Average coverage *3 cut off: " + str(arr_avg_cov[i-2]) + "\n")
            info_file.close()
            print("Sample " + str(i) +  "done.")

path = argv[1]
chr_name= argv[2]
pen_ref=int(argv[3])
pen_ind=int(argv[4])
arr_avg= argv[5]

ref_index= argv[6].split(',')
map_object = map(int,ref_index)
ref_index = list(map_object)

sample_index=(argv[7].split(','))
map_object = map(int,sample_index)
sample_index = list(map_object)
arr_avg_cov=calculate_avg_thr(arr_avg)

predict_deletions(path, chr_name, pen_ref, pen_ind, arr_avg_cov, ref_index, sample_index)
