# CumulativeClassifier
# Language: R
# Input: TXT (parameters)
# Output: CSV (optimal parameters)
# Tested with: PluMA 1.0, R 3.2.5

PluMA plugin that attempts to classify a set of samples 
using a specific set of features, merging as input both 
a clinical and a training set.  Both input files (training and clinical) are assumed to be in TSV (tab-separated value)
format, with rows corresponding to subjects and columns corresponding to feature values.
The input for the plugin is a TXT file of keyword-value pairs, training and clinical,
each mapping to their respective filenames.

This plugin is designed for time-series data and therefore also takes as input the selected 
features at each timepoint, with files in the format:
ChosenGenes_Cumulative_DEE1 RSV_(time).csv

where (time) is the current timepoint.  These consist of multiple lines of the form:

<featurename>,<value>

Optimal classification parameters for each timepoint are then sent to the output CSV
file provided to the plugin.  Timepoints are separated by empty string delimiters.

Note the input TSV and CSV files in the example/ directory are not publically available.
A future goal is to make a synthetic data set.  In the meantime however, one may
use it on their own tab-separated input data.
