# Make a python function which can compare output from script to the given .bed files and calculate:

# Accuracy = (TP+TN)(/TN+FN+TP+FP)
# Sensitivity (True Positive Rate) = TP/(TP+FN)
# Specificity (True Positive Rate) = TN/(FP+TN)
# ROC-curve

# True Positive - A true deletion (found through script and bed file)
# False Positive - Script find deletion but not in bed file
# True Negative - Icke deletion (interval from deletion to deletion)
# False Negative - Script does not find true deletion

# Find the smallest deletion (true positive with the smallest length)

def function1(bed1, bed2):
   # Compare 2 bed-files or two arrays or matrices
   # bed1 is from our script
   # bed2 is bed-file with true deletions 

   # True positives
   true_pos = 10
   # False positives
   false_pos = len(bed1) - true_pos
   # True negative
   true_neg = 12
   # False negative
   false_neg = len(bed2) - true_pos

   accuracy = (true_pos + true_neg)/(true_neg + false_neg + true_pos + false_pos)
   sensitivity = true_pos/(true_pos + false_neg)
   specificity = true_neg/(true_neg + false_pos)

   return accuracy, sensitivity, specificity

# Testing
bed1 = [30,36,80,98,195,200,345,455,560,600,602,690,740,765]
bed2 = [30,36,81,98,195,200,340,466,550,599,602,690,740,765,777,789,812,832]

print(function1(bed1,bed2))