import sys
import numpy as np
from hmmlearn import hmm
import matplotlib.pyplot as plt
import argparse
import pathlib

# To run script use following:
# zcat file.zip | head -n 100000 | python SimpleModel.py --sapmles 2 --fit 0 1 2 3 4 --plot

# Define the possible arguments
parser = argparse.ArgumentParser(description="Parse some arguments")
parser.add_argument('--samples', metavar='S', type=int, nargs='+', help='Samples to predict', required=False)
parser.add_argument('--samplestext', metavar='STex', type=pathlib.Path, nargs=1, help='Samples to predict', required=False)
parser.add_argument('--fit', metavar='F', type=int, nargs='+', help='Samples to fit model on', required=False)
parser.add_argument('--fittext', metavar='FTex', type=pathlib.Path, nargs=1, help='Samples to fit model on', required=False)
parser.add_argument('--plot', default=False, action='store_true')

# Hangle for parsed arguments
args = parser.parse_args()

# ---Assertion about the arguments---
# Samples to be predicted must be given by either the argument samples or samplestext
assert bool(args.samplestext) != bool(args.samples), "Both samples and samplestext can't be given."
# Samples to be fit thet model must be given by either the argument fit or fittext
assert not (bool(args.fittext) and bool(args.fit)), "Both fit and fittext can't be given."

#---Get data from pipe---
# Returns a 2D numpy array of the data
def get_data():
    data_pipeline = []
    for line in sys.stdin:
        # Append a list
        data_pipeline.append(line.split())
    data_pipeline = np.array(data_pipeline)
    return data_pipeline#.astype(np.int)

# Laod data from pipe
# X have shape: Position x Sample
sampleData = get_data()
chromosomeNames = sampleData[:, 0]
positions = sampleData[:, 1].astype(np.int)
reference = sampleData[:, 2].astype(np.int)
sampleData = sampleData[:, 3:].astype(np.int)
numOfPositions, numOfSamples = sampleData.shape

print(f"Number of samples: {numOfSamples}")

#---Samples to fit model to---
# Find the list of samples to predict
if args.fittext:
    with args.fittext[0].open() as file:
        # Read the file and then convert string to int
        samplesToFit = list(map(int, file.readlines()))
elif args.fit:
    samplesToFit = args.fit
else:
    # In case we train on all samples
    samplesToFit = list(range(numOfSamples))

#---Samples to predict---
# Find the list of samples to predict
if args.samplestext:
    with args.samplestext[0].open() as file:
        # Read the file and then convert string to int
        samplesToPredict = list(map(int, file.readlines()))
elif args.samples:
    samplesToPredict = args.samples

# Assert that the samples to predict is not out of bound by the data given
assert max(samplesToPredict) <= numOfSamples, "The sampels requested to predict is out of bound."


#---Define and fit model---
# Define HMM
model = hmm.GaussianHMM(n_components=2, covariance_type='tied', n_iter=100, tol=0.0001, verbose=True)
# Manipulate data for HMM, need when fitting multiple sequences
samplesFit = np.concatenate([sampleData[:,i] for i in samplesToFit])
lengths = [len(sampleData[:,i]) for i in samplesToFit]

model.fit(samplesFit.reshape(-1, 1), lengths)

#---Predic the states---
Z = []
for sample in samplesToPredict:
    Z.append(model.predict(sampleData[:, sample].reshape(-1, 1)))


#---Find the start of all states---
stateStarts = []
for i in range(len(samplesToPredict)):
    sampleStateStart = []
    for pos in range(numOfPositions-1):
        if Z[i][pos] != Z[i][pos+1]:
            sampleStateStart.append(pos)
    stateStarts.append(sampleStateStart[:]) # [:] is needed for deep copy

#---Calculate the average cover for each state---
stateAverageCover = []
for i in range(len(samplesToPredict)):
    sampleStateAverageCover = []
    pos = 0
    for stateStart in stateStarts[i]+[numOfPositions]:
        averageCoverage = 0
        posCounter = 0
        while pos < stateStart:
            averageCoverage += sampleData[pos, samplesToPredict[i]]
            pos += 1
            posCounter += 1
        averageCoverage = averageCoverage/posCounter
        sampleStateAverageCover.append(averageCoverage)
    stateAverageCover.append(sampleStateAverageCover[:])


# Implement cut-off at 20.
#for i in range(len(samplesToPredict)):


#print(stateStarts[0])
#print(stateAverageCover[0])

#plt.plot(stateAverageCover[0])
#plt.show()

# Plot the resulting predictions
if args.plot:
    for i in range(len(samplesToPredict)):
        plt.title(f"predicted state of sample: {samplesToPredict[i]}")
        plt.plot(Z[i])
        plt.show()
