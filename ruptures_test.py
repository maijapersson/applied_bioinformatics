#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 10:59:10 2021

@author: linnbeckman
"""
import pybedtools
import ruptures as rpt
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from remove_artifacts import roc_curve


#x100_covergae = pd.read_csv('/Users/linnbeckman/OneDrive/Uppsala Universitet_OneDrive/project_bioinf/Simulated_data/100X.coverage',header=None, sep='\t')


x100_cov_1mil = pd.read_csv('/Users/linnbeckman/OneDrive/Uppsala Universitet_OneDrive/project_bioinf/Simulated_data/1mil_100x', header=None, sep='\t')
x100_cov_300k = pd.read_csv('/Users/linnbeckman/OneDrive/Uppsala Universitet_OneDrive/project_bioinf/Simulated_data/300k_100x',header=None, sep='\t')


ind5=np.array(x100_cov_1mil[5])
ref=np.array(x100_cov_1mil[2])
true_bkps_bed = pybedtools.BedTool("1mil_C.bed")
#
# plt.figure(1)
# model = "l2"  # "l1", "rbf", "linear", "normal", "ar"
# algo = rpt.BottomUp(model=model, min_size=100).fit(ind5)
# ind5_bkps = algo.predict(pen=7000000)
# # rpt.show.display(ind5, ind5_bkps, figsize=(100, 3))
# # plt.title('Change Point Detection: Bottom Up Search Method')
# # plt.show()
# print(ind5_bkps)
#
# plt.figure(2)
# model = "l2"  # "l1", "rbf", "linear", "normal", "ar"
# algo = rpt.BottomUp(model=model, min_size=100).fit(ref)
# ref_bkps = algo.predict(pen=7000000)
# print(ref_bkps)
# # rpt.show.display(ref, ref_bkps, figsize=(100, 3))
# # plt.title('Change Point Detection: Bottom Up Search Method')
# # plt.show()
# print(ref_bkps)



#pred_bkps=pred_breakpoints(ref_bkps, ind5_bkps)
#print(pred_bkps)
#
# plt.figure(3)
# rpt.show.display(ind5, pred_bkps, figsize=(100, 3))
# plt.title('Change Point Detection: Bottom Up Search Method')
# plt.show()


ind5_bkps=[12320, 16475, 289420, 290640, 363280, 365350, 434810, 442620, 520990, 527215, 735955, 745230, 854485, 860590, 932730, 935905]
ref_bkps=[510735, 515010, 1000000]
roc_curve(100, ind5_bkps, true_bkps_bed, ind5, ref_bkps)
plt.show()



threshold= 5
avg_cov_arr=[]
pred_filterd=[]
# for i in range(len(ind5_bkps)-1):
#     start=ind5_bkps[i]
#     stop=ind5_bkps[i+1]
#
#     temp=ind5[start: stop]
#
#     avg_coverage=sum(temp)/len(temp)
#     if avg_coverage<threshold:
#         pred_filterd.extend([start, stop])
#
#     avg_cov_arr.append(avg_coverage)
#
#
# print(avg_cov_arr)
# print(pred_filterd)






# #Changepoint detection with the Pelt search method
# #20 000 took 10 min
# model="rbf"
# algo = rpt.Pelt(model=model, min_size=1000, jump=5).fit(points)
# result = algo.predict(pen=20)
# rpt.display(points, result, figsize=(10, 6))
# plt.title('Change Point Detection: Pelt Search Method')
# plt.show()
#
# #Changepoint detection with the Binary Segmentation search method
# model = "l2"
# algo = rpt.Binseg(model=model).fit(points)
# my_bkps = algo.predict(pen=10000)
# # show results
# rpt.show.display(points, my_bkps, figsize=(100, 6))
# plt.title('Change Point Detection: Binary Segmentation Search Method')
# plt.show()
#
#
# #Changepoint detection with window-based search method
# model = "l2"
# algo = rpt.Window(width=40, model=model).fit(points)
# my_bkps = algo.predict(n_bkps=10)
# rpt.show.display(points, my_bkps, figsize=(10, 6))
# plt.title('Change Point Detection: Window-Based Search Method')
# plt.show()
#
#
# #Changepoint detection with dynamic programming search method
# model = "l1"
# algo = rpt.Dynp(model=model, min_size=3, jump=5).fit(points)
# my_bkps = algo.predict(n_bkps=10)
# rpt.show.display(points, my_bkps, figsize=(10, 6))
# plt.title('Change Point Detection: Dynamic Programming Search Method')
# plt.show()
