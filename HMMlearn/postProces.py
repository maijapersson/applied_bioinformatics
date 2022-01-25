import json
import sys
import ruptures as rpt
from hmmlearn import hmm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pybedtools

# Takes output from hmm.py. Merges deletions that are close together and remove artifacts.

# Creates a bed object from the given array
def create_bed(array, chromosome):
    data_array = []
    if len(array) % 2 != 0:
        del array[-1]
    for i in range(0, len(array), 2):
        data_array.append([chromosome, array[i], array[i+1]])
    data_array = pd.DataFrame(data_array)
    return pybedtools.BedTool.from_dataframe(data_array)

# Remove all artifacts found in one reference from one sample.
def remove_artifacts_helper(reference, sample, chromosome):
    reference_bed = create_bed(reference, chromosome)
    sample_bed = create_bed(sample, chromosome)
    sample_bed = sample_bed.subtract(reference_bed, f=0.6,  A=True)
    return pybedtools.BedTool.to_dataframe(sample_bed).iloc[:, 1:3].to_numpy().flatten().tolist()

# Remove artifacts
# Makes changes in data_dict['sample_changes']
def remove_artifacts(data_dict):
    for i in range(len(data_dict['reference'])):
        for j in range(len(data_dict['sample'])):
            data_dict['sample_changes'][j] =remove_artifacts_helper(data_dict['reference_changes'][i], data_dict['sample_changes'][j], data_dict['chromosome'][j])
    return data_dict

# Plot the deeltions with coverage
def plot_deletions(data_dict):
    for i in range(len(data_dict['reference'])):
        rpt.display(np.array(data_dict['reference'][i]), np.array(data_dict['reference_changes'][i]))
        plt.title('Reference')
        plt.show()
    for i in range(len(data_dict['sample'])):
        rpt.display(np.array(data_dict['sample'][i]), np.array(data_dict['sample_changes'][i]))
        plt.title(f'sample {i}')
        plt.show()

# Merges deletions that are close to each other
# Makes changes in data_dict['reference_changes'] and data_dict['sample_changes']
# Looks at the stop position of a deletion and compaes it to the start of the next deletion
# Remember the datastructure data_dict[*_changes][i] have start and stop position for deletions and shape:
# [start, stop, start, stop, start, stop, ...] for a given reference/sample i
def merge_deletions(data_dict):
    for i in range(len(data_dict['reference_changes'])):
        # Ensure that deletions are present
        if data_dict['reference_changes'][i]:
            # List to hold the new deletions [start, stop, start ,stop, ...]
            new_changes = [data_dict['reference_changes'][i][0]]
            for pos in range(1, len(data_dict['reference_changes'][i])-1, 2):
                if (data_dict['reference_changes'][i][pos+1] - data_dict['reference_changes'][i][pos]) > 5000:
                    new_changes.append(data_dict['reference_changes'][i][pos])
                    new_changes.append(data_dict['reference_changes'][i][pos+1])
            new_changes.append(data_dict['reference_changes'][i][-1])
            data_dict['reference_changes'][i] = new_changes

    for i in range(len(data_dict['sample_changes'])):
        # Ensure that deletions are present
        if data_dict['sample_changes'][i]:
            # List to hold the new deletions [start, stop, start ,stop, ...]
            new_changes = [data_dict['sample_changes'][i][0]]
            for pos in range(1, len(data_dict['sample_changes'][i])-1, 2):
                if (data_dict['sample_changes'][i][pos+1] - data_dict['sample_changes'][i][pos]) > 5000:
                    new_changes.append(data_dict['sample_changes'][i][pos])
                    new_changes.append(data_dict['sample_changes'][i][pos+1])
            new_changes.append(data_dict['sample_changes'][i][-1])
            data_dict['sample_changes'][i] = new_changes

    return data_dict

# Save the reference and sample deletions as .bed files.
def save_bed_file(data_dict):
    for i in range(len(data_dict['reference'])):
        create_bed(data_dict['reference_changes'][i], data_dict['chromosome'][0]).saveas(f'reference_{i}.bed')
    for i in range(len(data_dict['sample'])):
        create_bed(data_dict['sample_changes'][i], data_dict['chromosome'][0]).saveas(f'sample_{i}.bed')


if __name__ == '__main__':
    for line in sys.stdin:
        data_dict = json.loads(line)

    data_dict = merge_deletions(data_dict)
    data_dict = remove_artifacts(data_dict)
#    plot_deletions(data_dict)
    save_bed_file(data_dict)
