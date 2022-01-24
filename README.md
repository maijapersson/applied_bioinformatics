# Project code in course Applied Bioinformatics

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
