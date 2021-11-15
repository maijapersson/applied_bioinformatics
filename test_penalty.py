import ruptures as rpt
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from remove_artifacts import pred_breakpoints
from benchmarking import function1
import pybedtools


x100_cov_300k = pd.read_csv('/Users/linnbeckman/OneDrive/Uppsala Universitet_OneDrive/project_bioinf/Simulated_data/1mil_100x',header=None, sep='\t')

true_del_bed = pybedtools.BedTool("300k_C.bed")
sample5=np.array(x100_cov_300k[5])
reference=np.array(x100_cov_300k[2])

def roc_curve(sample, ref, true):

    penalty=np.linspace(50000,5000000, 10)

    sensitivity=[]
    specificity=[]

    model = "l2"  # "l1", "rbf", "linear", "normal", "ar"
    for pen in penalty:
        algo = rpt.BottomUp(model=model, min_size=100).fit(sample)
        sample_bkps = algo.predict(pen=pen)

        algo = rpt.BottomUp(model=model, min_size=100).fit(ref)
        ref_bkps = algo.predict(pen=pen)

        pred_breakpoints(ref_bkps, sample_bkps)

        benchmarks=function1("removed_artifacts_pred.bed",true, 1000000)

        sensitivity.append(benchmarks[1])
        specificity.append(benchmarks[2])

    return sensitivity, specificity



#test

sens, spec = roc_curve(sample5, reference, true_del_bed)

fpr=1-np.array(spec)

print(sens)
print(spec)

plt.scatter(fpr, sens)
plt.show()
