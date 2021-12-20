#Display multiple samples

from sys import argv
import pandas as pd
import pybedtools
import matplotlib.pyplot as plt
import ruptures as rpt

import numpy as np
import pandas as pd


def bed_to_array(bed):
    return pybedtools.BedTool.to_dataframe(bed).loc[:, ["start", "end"]].to_numpy().flatten()

def display_changepoints(sample_index_arr, reference_index_arr, path, avg_arr):
    data=pd.read_csv(path, header=None, sep='\t')
    avg_cutoff=pd.read_csv(avg_arr, header=None)

    fig = plt.figure(num = "ruptures_figure"  , figsize=(100, 50))
    #define grid of nrows x ncols
    gs = fig.add_gridspec(len(sample_index_arr)+1, 1)
    start=0

    individual=np.array(data[2])
    individual[individual > avg_cutoff.at[0,0]*3] =avg_cutoff.at[0,0]*3
    bkps_bed=pybedtools.BedTool('results/reference.bed')
    ind_bkps=bed_to_array(bkps_bed)
    _, curr_ax = rpt.display(individual, ind_bkps, num="ruptures_figure")
    curr_ax[0].set_position(gs[start].get_position(fig))
    curr_ax[0].set_subplotspec(gs[start])
    plt.title('Reference after filtration')
    plt.xlabel('Position')
    plt.ylabel('Coverage')
    plt.xticks(np.arange(0, len(individual)+1, 100000))
    start=1

    for index in sample_index_arr:
        individual=np.array(data[index])
        individual[individual > avg_cutoff.at[index-2,0]*3] =avg_cutoff.at[index-2,0]*3
        bkps_bed=pybedtools.BedTool('results/sample_'+str(index)+'.bed')
        ind_bkps=bed_to_array(bkps_bed)

        _, curr_ax = rpt.display(individual, ind_bkps, num="ruptures_figure")
        curr_ax[0].set_position(gs[start].get_position(fig))
        curr_ax[0].set_subplotspec(gs[start])
        plt.title('Sample '+str(index)+'. After filtration and removal of artifacts')
        plt.xlabel('Position')
        plt.ylabel('Coverage')
        plt.xticks(np.arange(0, len(individual)+1, 100000))
        start=start+1
    fig.savefig('results/sample_figure_test.png')
    plt.show()
all=[3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20, 21]
display_changepoints(all, [2], '2mil', 'avg_arr.csv')
