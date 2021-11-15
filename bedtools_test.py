import ruptures as rpt
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from remove_artifacts import pred_breakpoints
from benchmarking import function1
import pybedtools

x100_cov_300k = pd.read_csv('/Users/linnbeckman/OneDrive/Uppsala Universitet_OneDrive/project_bioinf/Simulated_data/300k_100x',header=None, sep='\t')
true_bkps_bed = pybedtools.BedTool("300k_C.bed")
ind5=np.array(x100_cov_300k[5])
ref=np.array(x100_cov_300k[2])


def test_diff_pen_bottomUp(costfunction, ref, ind, true):
    #model = "l2"  # "l1", "rbf", "linear", "normal", "ar"

    pen=np.linspace(1000,1000000, 10)
    sensitivity=[]
    specificity=[]
    p=1000
    # for p in pen:
        #Fit the individual
    algo = rpt.BottomUp(model=costfunction, min_size=100).fit(ind)
    ind_bkps = algo.predict(pen=p)

        #Fit the reference
    algo = rpt.BottomUp(model=costfunction, min_size=100).fit(ref)
    ref_bkps = algo.predict(pen=p)

        #Remove artefacts
    pred_breakpoints(ref_bkps, ind_bkps)

        #Calculate benchmark values
    benchmarks=function1("removed_artifacts_pred.bed",true, 300000)

    sensitivity.append(benchmarks[1])
    specificity.append(benchmarks[2])



    return sensitivity, specificity


sens, spec =test_diff_pen_bottomUp("l2", ind5, ref, true_bkps_bed)

print(sens)
print(spec)

plt.plot(100-spec, sens)
