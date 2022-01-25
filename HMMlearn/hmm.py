import json
import sys
import ruptures as rpt
from hmmlearn import hmm
import numpy as np
import matplotlib.pyplot as plt
import time

# Reads the output from removeOutliers.py. Creates HMMs to find change points in
# data_dict["reference"] and data_dict["sample"] and classifies the found
# sequences as deletions or non-deletions based on their average coverage
# Adds the following keys:
# data_dict["reference_changes"] -> 2D list, for each reference contains a list with the change points of that raference. Ex [0, 100, 101, 142, 143, 210, ...]
# data_dict["sample"] -> 2D list, for each sample contains a list with the change points of that saample.
# data_dict['reference_del_avg_cov'] -> 2D list, for each reference contains a list with the average coverage for each deletion in that reference
# data_dict['sample_del_avg_cov'] -> 2D list, for each sample contains a list with the average coverage for each deletion in that sample


# Takes the output from a HMM and returns the positions where states start and stop
# returns list [start1, stop1, start2, stop2, ..., startn, stopn] if there are n changes
def convert_state(hmm_output):
    hmm_output[0] = hmm_output[1]
    hmm_output[-1] = hmm_output[-2]
    changes = [0]
    for i in range(len(hmm_output)-1):
        if hmm_output[i] != hmm_output[i+1]:
            changes.append(i)   # stop of this state
            changes.append(i+1) # start of next state
    changes.append(len(hmm_output)-1)
    if changes[1] == changes[0]:
        changes.pop(0)
    if changes[-1] == changes[-2]:
        changes.pop()
    return changes

# Calculate the average coverage of area from start to stop.
def calc_average_coverage(start, stop, coverage):
    return sum(coverage[start: stop+1]) / (stop+1-start)

# Trains HMMs on the reference and sample coverage
def train_models(data_dict):
    # Create the HMMs
    reference_hmm = [hmm.GaussianHMM(n_components=2, n_iter=100, covariance_type='tied', tol=1, verbose=False)] * len(data_dict['reference'])
    sample_hmm = [hmm.GaussianHMM(n_components=2, n_iter=100, covariance_type='tied', init_params='' ,tol=1, verbose=False)] * len(data_dict['sample'])
    # Train the HMMs
    for i in range(len(data_dict['reference'])):
        # Have not tested of the min() function workds out, ie. if a sample of more than 10000000 bp is run it will only train on the first 10000000.
        reference_hmm[i].fit(np.array(data_dict['reference'][i][0:min(10000000, len(data_dict['reference'][i]))]).reshape(-1, 1))
    for i in range(len(data_dict['sample'])):
        sample_hmm[i].fit(np.array(data_dict['sample'][i][0:min(10000000, len(data_dict['sample'][i]))]).reshape(-1, 1))
    # Add the HMM to the dictionary, must be removed before the dictionary is printed with json.
    data_dict['reference_hmm'] = reference_hmm
    data_dict['sample_hmm'] = sample_hmm
    return data_dict

# Uses the trained HMMs to predict the state of each position
# Creates the keys data_dict["reference_changes"] and data_dict["sample_changes"]
def predict_state(data_dict):
    # Initiate list of state changes
    reference_deletion = [0] * len(data_dict['reference'])
    sample_deletion = [0] * len(data_dict['sample'])
    # Predict the references
    for i in range(len(data_dict['reference'])):
        predicted_states = data_dict['reference_hmm'][i].predict(np.array(data_dict['reference'][i]).reshape(-1, 1))
        reference_deletion[i] = convert_state(predicted_states)
    data_dict['reference_changes'] = reference_deletion
    # Predict the samples
    for i in range(len(data_dict['sample'])):
        predicted_states = data_dict['sample_hmm'][i].predict(np.array(data_dict['sample'][i]).reshape(-1, 1))
        sample_deletion[i] = convert_state(predicted_states)
    data_dict['sample_changes'] = sample_deletion
    return data_dict

# Applies the cutoff filter and classifies deletions
# Removes the areas in data_dict['reference_changes'] and data_dict['sample_changes'] that have an average coverage too high
# Creates the keys data_dict['reference_del_avg_cov'] and data_dict['sample_del_avg_cov']
def cutoff_filter(data_dict):
    data_dict['reference_del_avg_cov'] = []
    data_dict['sample_del_avg_cov'] = []
    # For reference remove areas and calculate average coverage for remaining deletions
    for i in range(len(data_dict['reference'])):
        average_coverage = []
        deletion_pos = [] # [start, stop]
        for j in range(len(data_dict['reference_changes'][i])//2):
            start = data_dict['reference_changes'][i][j*2]
            stop = data_dict['reference_changes'][i][j*2+1]
            avg_cov = calc_average_coverage(start, stop, data_dict['reference'][i])
            if avg_cov < data_dict['reference_cutoff'][i]:
                average_coverage.append(avg_cov)
                deletion_pos.append(start)
                deletion_pos.append(stop)
        data_dict['reference_changes'][i] = deletion_pos
        data_dict['reference_del_avg_cov'].append(average_coverage)

    # For sample remove areas and calculate average coverage for remaining deletions
    for i in range(len(data_dict['sample'])):
        average_coverage = []
        deletion_pos = [] # [start, stop]
        for j in range(len(data_dict['sample_changes'][i])//2):
            start = data_dict['sample_changes'][i][j*2]
            stop = data_dict['sample_changes'][i][j*2+1]
            avg_cov = calc_average_coverage(start, stop, data_dict['sample'][i])
            if avg_cov < data_dict['sample_cutoff'][i]:
                average_coverage.append(avg_cov)
                deletion_pos.append(start)
                deletion_pos.append(stop)
        data_dict['sample_changes'][i] = deletion_pos
        data_dict['sample_del_avg_cov'].append(average_coverage)
    return data_dict

# Removes the HMMs from data_dict, erquired before printing data with json
def remove_hmm(data_dict):
    del data_dict['reference_hmm']
    del data_dict['sample_hmm']
    return data_dict

# Plot stuff
def plot_deletions(data_dict):
    for i in range(len(data_dict['reference'])):
        rpt.display(np.array(data_dict['reference'][i]), np.array(data_dict['reference_changes'][i]))
        plt.title('Reference')
        plt.show()
    for i in range(len(data_dict['sample'])):
        rpt.display(np.array(data_dict['sample'][i]), np.array(data_dict['sample_changes'][i]))
        plt.title(f'sample 3')
        plt.show()


if __name__ == '__main__':
    for line in sys.stdin:
        data_dict = json.loads(line)

    # Trains HMM models, these are saved in data_dict
    data_dict = train_models(data_dict)
    # Uses those HMM models to classify changes
    data_dict = predict_state(data_dict)
    # Classify the areas between changes as either deletion or non-deleteion
    data_dict = cutoff_filter(data_dict)
    # Remove the HMMs since they can not be passed on by json
    data_dict = remove_hmm(data_dict)
    print(json.dumps(data_dict))
