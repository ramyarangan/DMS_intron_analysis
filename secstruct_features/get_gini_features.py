"""
Compute Gini index in windows for a sequence
"""

from os import path
import numpy as np
import os
from matplotlib import pyplot as plt 
from scipy import stats 
import pandas as pd 
import seaborn as sns 
from secstruct_util import *

MUTFREQ_COVERAGE = 6456 # 1971 * 2

def gini(reacs):
	# Taken from https://neuroplausible.com/gini
	"""Calculate the Gini coefficient of a numpy array."""
	# All values are treated equally, arrays must be 1d:
	reacs = reacs.flatten()
	# Values cannot be 0:
	reacs += 0.0000001
	# Values must be sorted:
	reacs = np.sort(reacs)
	# Index per array element:
	index = np.arange(1,reacs.shape[0]+1)
	# Number of array elements:
	n = reacs.shape[0]
	# Gini coefficient:
	return ((np.sum((2 * index - n  - 1) * reacs)) / (n * np.sum(reacs)))

def get_gini_windows(name, seq, mut_freq_file, window_size=20, window_scan=10, \
	overhang_5p=0, overhang_3p=0):
	gini_vals = []
	num_passing = 0

	if name.replace('\n', '') in snorna_introns:
		return []

	scan_start = overhang_5p
	scan_end = len(seq) - overhang_3p

	freqs, _, totals, seq = get_mut_freq(mut_freq_file, name)

	freqs = freqs[scan_start:scan_end]
	totals = totals[scan_start:scan_end]
	seq = seq[scan_start:scan_end]

	seq = np.array(list(seq))
	seq_mask = np.logical_or(seq == "A", seq == "C")
	freqs = freqs[seq_mask]
	totals = totals[seq_mask]

	start_idxs = np.arange(0, len(freqs) - window_size, window_scan)
	for start_idx in start_idxs:
		cur_reac_window = freqs[start_idx:(start_idx + window_size)]
		cur_min_totals = min(totals[start_idx:(start_idx + window_size)])
		if cur_min_totals > MUTFREQ_COVERAGE:
			gini_vals += [gini(cur_reac_window)]

	return gini_vals

def get_gini_windows_tag(tag, fasta_file, mut_freq_file, window_size=20, window_scan=10):
	names, seqs, name_seq_dict, _ = get_names_seqs_from_fasta(fasta_file)
	gini_vals = []
	num_passing = 0
	for name in names: 
		if name != tag:
			continue
		if name.replace('\n', '') in snorna_introns:
			return gini_vals

		return get_gini_windows(name, name_seq_dict[name], mut_freq_file, \
			window_size=window_size, window_scan=window_scan)

	return gini_vals
