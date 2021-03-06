# Project code in course Applied Bioinformatics
##### Authors: Linn Beckman, Maija Persson and Viktor Öhrner

## How to run finding_deletions_ruptures.py

Module dependencies: ruptures, pybedtools, pandas, numpy, matplotlib.pyplot, time, sys and argv.

The script requires that a dictionary with the name results is created at the same place as the file finding_deletions_ruptures.py.

**The script can run from the command line by the following command:**

> python finding_deletions_ruptures.py path_to_detection_file chr_name pen_ref pen_ind path_to_filtration_file ref_index sample_index

**EXAMPLE:**
> python finding_deletions_ruptures.py 10Mb_coverage.csv Chr01 100000 100000 30Mb_coverage.csv 2,3 4,5,6

*path_to_detection_file:* The path to the file which contains coverage information for all individuals including references which the change point detection analysis should be performed on. 

*chr_name:* The chromosome name which will be shown in the resulting bed files.

*pen_ref:* Penalty used for the references, recommended: 100000.

*pen_ind:* Penalty used for the samples,  recommended: 100000.

*path_to_filtration_file:* The path to the coverage file which is used for calculating average coverage in the samples and deletion threshold (all samples and references). OBS! needs a file with at least 30 Mb.

*ref_index:* Indexes for the reference samples in the coverage file. 

*sample_index (optional):* Indexes for the samples which you want to analyze. If not specified, then all samples except the specified reference(s) will be analyzed.

**OUTPUTS:**

All outputs are stored in the results dictionary. The script outputs one bed file with the start and stop positions of the detected deletion/artefacts for each sample or reference. For each sample and reference, a plot with the resulting change points are saved. There is also a text file for each sample and reference containing additional information about the analysis and a plot with all samples and references.


## How to run the HMM pipeline

Module dependencies: hmm_learn, scikit-learn, numpy, matplotlib.pyplot, scipy, json, pathlib, ruptures, pybedtools, pandas, collections, sys and argparse, time, pathlib.

**The script can run from the command line by the following command:**

> zcat DATA_PATH/FILE.gz | head -n INT | python createDict.py --samples FILE1.tex --references FILE2.tex | python removeOutliers.py | python hmm.py | python postProces.py

Where FILE1.tex and FILE2.tex have the index of the samples and references on each row. These .py files can be found in the folder ./HMMlearn with examples of FILE1.tex and FILE2.tex.

**EXAMPLE:**

> zcat DATA_PATH/PA_chr01.depth.gz | head -n 1000000 | python createDict.py --sample sample.tex --reference reference.tex | python removeOutliers.py | python hmm.py | python postProces.py

**OUTPUTs:**

Outputs a .bed file with deletions for each reference and sample. A plot with deletions and coverage can be generated as well. 
