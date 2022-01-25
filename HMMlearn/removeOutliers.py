import json
import sys
import collections

# This script reads the output from createDict.py. Brings down the values in
# data_dict['reference'] and data_dict["sample"] that exceed the average
# coverage * 3.
# Adds the following keys to the dictionary:
# data_dict['reference_average'] -> list of average coverage of all references
# data_dict['sample_average'] -> list of average coverage of all samples
# data_dict['reference_cutoff'] -> list of the cut off value for the references
# data_dict['sample_cutoff'] -> list of the cut off value for the samples

# Takes array of coverage
def calc_average(coverage):
    occurance = collections.Counter(coverage)
    # Delete coverage 0-4
    for i in range(5):
        del occurance[i]
    # Delete the last coverage
    del occurance[max(occurance.keys())]
    return max(occurance, key=occurance.get)

def removeOutliers(data_dict):
    num_pos = len(data_dict['reference'][0])
    data_dict['reference_average'] = [0] * len(data_dict['reference'])
    data_dict['sample_average'] = [0] * len(data_dict['sample'])

    for i in range(len(data_dict['reference'])):
        data_dict['reference_average'][i] = calc_average(data_dict['reference'][i])
    for i in range(len(data_dict['sample'])):
        data_dict['sample_average'][i] = calc_average(data_dict['sample'][i])

    # Calculate cut-off
    reference_cutoff = [data_dict['reference_average'][i] * 3 for i in range(len(data_dict['reference']))]
    sample_cutoff = [data_dict['sample_average'][i] * 3 for i in range(len(data_dict['sample']))]

    # Applies cut off
    #print("Start appyling reference cut-off")
    for i in range(len(data_dict['reference'])):
        data_dict['reference'][i] = [min(float(data_dict['reference'][i][pos]), reference_cutoff[i]) for pos in range(len(data_dict['reference'][i]))]

    #print("Start appyling sample cut-off")
    for i in range(len(data_dict['sample'])):
        data_dict['sample'][i] = [min(float(data_dict['sample'][i][pos]), sample_cutoff[i]) for pos in range(len(data_dict['sample'][i]))]

    return data_dict


def find_cutoff(data_dict):
    reference_cutoff = [0] * len(data_dict['reference'])
    sample_cutoff = [0] * len(data_dict['sample'])

    for i in range(len(data_dict['reference'])):
        reference_cutoff[i] = data_dict['reference_average'][i] * 0.55
    for i in range(len(data_dict['sample'])):
        sample_cutoff[i] = data_dict['sample_average'][i] * 0.15

    data_dict['reference_cutoff'] = reference_cutoff
    data_dict['sample_cutoff'] = sample_cutoff

    return data_dict

if __name__ == '__main__':
    for line in sys.stdin:
        data_dict = json.loads(line)

    # Remove outliers
    data_dict = removeOutliers(data_dict)
    data_dict = find_cutoff(data_dict)

    print(json.dumps(data_dict))
