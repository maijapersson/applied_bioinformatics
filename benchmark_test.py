#!/usr/bin/env python3
# -*- coding: utf-8 -*-


#producing some plots that show the % of missed sites,
#% of over-estimated sites,
#and also looking at deletions of different length and how this performs

import pybedtools
import pandas as pd
import sys
import matplotlib.pyplot as plt
import numpy as np

pred_bkps_bed = pybedtools.BedTool("removed_artifacts_pred.bed")
true_bkps_bed = pybedtools.BedTool("1mil_C.bed")

true_del_with_matches = true_bkps_bed.intersect(pred_bkps_bed, c=True).to_dataframe()
true_pos= pred_bkps_bed.intersect(true_bkps_bed)

true_pos=pred_bkps_bed.intersect(true_bkps_bed).to_dataframe()
true_del=true_del_with_matches

print("----True del ------")
print(true_del)
print("----True positives  ------")
print(true_pos)

false_pos=pred_bkps_bed.subtract(true_bkps_bed).to_dataframe(header=None)
print("----- False pos ------")
print(false_pos)
# print("-------True pos----------")
# print(true_pos)
# print("-------True deletion----------")
# print(true_bkps_bed)
# print("---------Pred deletion--------")
# print(pred_bkps_bed )
#
# ms_intersect=true_pos.intersect(true_bkps_bed, F=1)#Exact matches
# missed_sites= true_bkps_bed.subtract(true_pos)# False negatives
#
# print("MS_INTERSECT")
# print(ms_intersect)


# feature = true_pos[0]
# true_pos_arr=[]
# false_pos_arr=[]
#
# missed_sites_arr=[]
#
# for feature in missed_sites:
#     missed_sites_arr.append(feature.stop-feature.start)
# print("Missed missed_sites")
# print(missed_sites)
#
# #
# # for feature in true_pos:
# #     true_pos_arr.append(feature.stop-feature.start)
# #
# # for feature in false_pos:
# #     false_pos_arr.append(feature.stop-feature.start)
# #
real_del=[]
for feature in true_bkps_bed:
    real_del.append(feature.stop-feature.start)

missed_sites = []
for i in range(true_del.shape[0]):
    if true_del.loc[i,"name"]==0:
        missed_sites.append(true_del.loc[i,"end"]-true_del.loc[i,"start"])
    else:
        j = 0
        while true_pos.loc[j,"start"] < true_del.loc[i,"end"]:
            j=j+1
            if j >= true_pos.shape[0]:
                break
        missed_sites.append(true_pos.loc[j-1,"start"]-true_del.loc[i,"start"] + true_del.loc[i,"end"]-true_pos.loc[j-1,"end"])
# test_missed_sites=[]
#
# for i in range(true_del.shape[0]):
#     start = true_pos.loc[i,"start"]
#     stop = true_pos.loc[i,"end"]
#     if true_del.loc[i,"name"]==0:
#         test_missed_sites.append(true_del.loc[i,"end"]-true_del.loc[i,"start"])
#     else:
#         test_missed_sites.append(start-true_del.loc[i,"start"] + true_del.loc[i,"end"]-stop)



print("-------missde sites --------")
print(missed_sites)
# division of lists
# using zip() + list comprehension

res = [i / j for i, j in zip(missed_sites, real_del)]
print("----------real deletions ----------")
print(real_del)
# printing result
print ("The division list is : " + str(res))

x=real_del

plt.scatter(x, res)


print(res)
#

#
# true_neg = 300000 - sum(false_pos_arr) - false_neg -sum(true_pos_arr)#real non del - false pos = (total lenght - real del) - false pos
# # False negative
# #Real del - true pos
#
# print("True positive: "+ str(sum(true_pos_arr)))
# print("False positive: " + str(sum(false_pos_arr)))
# print("False negative: " + str(false_neg))
# print("True negative: " + str(true_neg))

# Script for calculating over estimation of sites for true deletions
# Skips false positives with no correlation to true deletions
over_est_sites = []
# Loop through size/number of true deletions
for i in range(true_del.shape[0]):
    left_over_est = 0 # Over estimation to the left of the deletion
    right_over_est = 0 # Over estimation to the right of the deletion
    # Loop through all false positives
    for j in range(false_pos.shape[0]):
        # Left side match
        if false_pos.loc[j,"end"] == true_del.loc[i,"start"]:
            left_over_est = false_pos.loc[j,"end"]-false_pos.loc[j,"start"]
        # Right side match
        elif false_pos.loc[j,"start"] == true_del.loc[i,"end"]:
            right_over_est = false_pos.loc[j,"end"]-false_pos.loc[j,"start"]
    over_est_sites.append(left_over_est + right_over_est) # Append number of over estimated sites

print("Over estimates")
print(over_est_sites)

p_over_est = [i / j for i, j in zip(over_est_sites, real_del)]

print(p_over_est)
plt.scatter(real_del, p_over_est)
plt.show()

# over_estimated_sites=[] #False positives
#
# for fp in range(false_pos.shape[0]):
#     for td in range(true_del.shape[0]):
#         if false_pos.loc[fp,"start"]==true_del.loc[td,"end"] or false_pos.loc[fp,"end"]==true_del.loc[td,"start"]:
#             j=j+1
#             over_estimated_sites.append(false_pos.loc[fp,"end"]-false_pos.loc[fp,"start"])
#     i=i+1
