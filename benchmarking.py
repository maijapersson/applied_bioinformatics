# Make a python function which can compare output from script to the given .bed files and calculate:

# Accuracy = (TP+TN)(/TN+FN+TP+FP)
# Sensitivity (True Positive Rate) = TP/(TP+FN)
# Specificity (True Positive Rate) = TN/(FP+TN)
# ROC-curve

# True Positive - A true deletion (found through script and bed file)
# False Positive - Script find deletion but not in bed file (over estimated)
# True Negative - Icke deletion (interval from deletion to deletion)
# False Negative - Script does not find true deletion (missed sites)

# Find the smallest deletion (true positive with the smallest length)

import pybedtools

def calc_acc_sens_spec(pred_bkps_bed, true_bkps_bed):
   # Compare 2 bed-files or two arrays or matrices
   # pred_bed is from our script
   # true_bed is bed-file with true deletions

   # pred_bkps_bed = pybedtools.BedTool(pred_bed)
   # true_bkps_bed = pybedtools.BedTool(true_bed)

   true_pos_bed= pred_bkps_bed.intersect(true_bkps_bed)
   false_pos_bed=pred_bkps_bed.subtract(true_bkps_bed)

   true_pos_arr=[]
   false_pos_arr=[]
   real_del=[]

   for feature in true_pos_bed:
       true_pos_arr.append(feature.stop-feature.start)
   for feature in false_pos_bed:
       false_pos_arr.append(feature.stop-feature.start)
   for feature in true_bkps_bed:
       real_del.append(feature.stop-feature.start)

   # True positives
   true_pos = sum(true_pos_arr)

   # False positives
   false_pos =  sum(false_pos_arr)
   # False negative
   false_neg = sum(real_del)- sum(true_pos_arr)#Real deletions - true positives
   # True negative
   true_neg = 1000000 - sum(false_pos_arr) - false_neg - sum(true_pos_arr)


   accuracy = (true_pos + true_neg)/(true_neg + false_neg + true_pos + false_pos)
   sensitivity = true_pos/(true_pos + false_neg)
   specificity = true_neg/(true_neg + false_pos)

   return accuracy, sensitivity, specificity

def over_estimated_sites(true_del, false_pos):
    over_est_sites = []
    real_del=[]
    # Loop through size/number of true deletions
    for i in range(true_del.shape[0]):
        real_del.append(true_del.loc[j,"end"]-true_del.loc[j,"start"])
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

    p_over_est = [i / j for i, j in zip(over_est_sites, real_del)]

    return p_over_est



def under_estimated_sites(true_del, true_pos):
    missed_sites = []
    real_del=[]
    for i in range(true_del.shape[0]):
        real_del.append(true_del.loc[j,"end"]-true_del.loc[j,"start"])
        if true_del.loc[i,"name"]==0:
            missed_sites.append(true_del.loc[i,"end"]-true_del.loc[i,"start"])
        else:
            j = 0
            while true_pos.loc[j,"start"] < true_del.loc[i,"end"]:
                j=j+1
                if j >= true_pos.shape[0]:
                    break
            missed_sites.append(true_pos.loc[j-1,"start"]-true_del.loc[i,"start"] + true_del.loc[i,"end"]-true_pos.loc[j-1,"end"])

    fractions = [i / j for i, j in zip(missed_sites, real_del)]
    return fractions
