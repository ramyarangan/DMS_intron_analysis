"""
Evaluate DMS reactivity for rRNA residues:
1. Get / plot correlation with Zubradt, et al. (2017) rRNA reactivity values
2. Get ROC curve for classifying rRNA base-paired vs unpaired + accessible residues
"""
from sklearn import metrics
from scipy.stats import pearsonr
import numpy as np
from matplotlib import pyplot as plt 
import xml.etree.ElementTree as ET
import math 
from os import path

EXP = 'd45'

secstruct_filename_dict = {
	'18s': 'rrna_secstruct/18s_4v88_secstruct.txt',
	'25s': 'rrna_secstruct/25s_4v88_secstruct.txt'
}

fasta_name_dict = {
	'18s': 'RDN18-1',
	'25s': 'RDN25-1'
}

coverage_cutoffs = {
	'd45': 1097 # Read depth required for 20 RPKM 
}

rrna_len_dict = {
	'18s': 1800
}

# Read in pre-computed Solvent Accessible Surface Area data 
# (computed in PyMOL in with script rrna_secstruct/)
def get_is_accessible(sasa_filename, sasa_cutoff=2):
	f = open(sasa_filename)
	sasa_lines = f.readlines()
	f.close()

	is_accessible = []
	for sasa_line in sasa_lines:
		sasa_val = float(sasa_line.split()[1])
		is_accessible += [sasa_val > sasa_cutoff]

	return is_accessible

# Reads in secondary structure file (obtained from DSSR on xtal structure)
# Returns a list of -1 if no info, 0 if not base-paired and accessible, 1 if base-paired
def get_is_base_paired_xtal(rrna_name):
	is_accessible = []
	sasa_filename = 'rrna_secstruct/sasa_dotdens1_solvradius3_aN1_cN3_bundle12_' + rrna_name + '.txt'
	is_accessible = get_is_accessible(sasa_filename) # sasa_dotdens1_solvradius3_aN1_cN3_18s
	
	secstruct_filename = secstruct_filename_dict[rrna_name]
	f = open(secstruct_filename)
	lines = f.readlines()
	f.close()

	seq = lines[0].replace('\n', '')
	all_chars = list(lines[1].replace('\n', ''))

	is_base_paired = []
	for ii, cur_char in enumerate(all_chars):
		if not ((seq[ii] == 'a') or (seq[ii] == 'c')):
			is_base_paired += [-1]
		elif cur_char == '-':
			is_base_paired += [-1]
		elif cur_char == '.':
			if is_accessible[ii]:
				is_base_paired += [0]
			else:
				is_base_paired += [-1]
		else:
			# captures (, ), {, }, [, ]
			is_base_paired += [1]
	return np.array(is_base_paired)

# Read in mutational frequencies from DMS-MaPseq
def get_mut_freq(mut_freq_file, fasta_name):
	f = open(mut_freq_file)
	mut_freq_lines = f.readlines()
	f.close()

	freqs = []
	muts = []
	totals = []
	for ii in range(int(len(mut_freq_lines)/5)):
		name = mut_freq_lines[ii*5].replace('\n', '')
		if name == fasta_name:
			seq = mut_freq_lines[ii*5 + 1].replace('\n', '')
			seq = np.array(list(seq))

			muts = mut_freq_lines[ii*5 + 2].split(',')
			muts = np.array([int(x) for x in muts])

			totals = mut_freq_lines[ii*5 + 3].split(',')
			totals = np.array([int(x) for x in totals])

			freqs = muts/totals
	return np.array(freqs), np.array(muts), np.array(totals), seq

# Get base-paired positions based on mutational frequencies with a given
# cutoff for classifying paired from unpaired
def get_is_base_paired_mutfrac(mut_freq_file, fasta_name, cur_cutoff):
	f = open(mut_freq_file)
	mut_freq_lines = f.readlines()
	f.close()

	is_reactive = []

	for ii in range(int(len(mut_freq_lines)/5)):
		name = mut_freq_lines[ii*5].replace('\n', '')
		if name == fasta_name:
			muts = mut_freq_lines[ii*5 + 2].split(',')
			muts = np.array([int(x) for x in muts])

			totals = mut_freq_lines[ii*5 + 3].split(',')
			totals = np.array([int(x) for x in totals])

			freqs = muts/totals

			is_reactive = freqs > cur_cutoff

	is_base_paired = []
	for cur_base in is_reactive:
		if not cur_base:
			is_base_paired += [1]
		else:
			is_base_paired += [0]
	return np.array(is_base_paired)

# Get correlation for mutational frequencies: our data vs Zubradt, et al. (2017)
def plot_rrna_correlation(rrna_names = ['18s', '25s']):
	mut_freq_file = "combined_1221/rfcount/" + EXP + "_rrna_all_view.txt"
	rouskin_mut_freq_file = "combined_1221/rfcount/rouskin_rrna_view.txt"

	all_freqs = []
	all_rouskin_freqs = []
	for rrna_name in rrna_names:
		fasta_name = fasta_name_dict[rrna_name]
		freq, muts, totals, seq = get_mut_freq(mut_freq_file, fasta_name)
		rouskin_freq, rouskin_muts, rouskin_totals, seq = \
			get_mut_freq(rouskin_mut_freq_file, fasta_name)

		coverage_cutoff = coverage_cutoffs[EXP]
		coverage_mask = np.logical_and(totals > coverage_cutoff, \
			rouskin_totals > coverage_cutoff)
		seq_mask = np.logical_or((seq == 'A'), (seq == 'C')) 
		mask = np.logical_and(seq_mask, coverage_mask)
		freq = freq[mask]
		muts = muts[mask]
		totals = totals[mask]
		rouskin_freq = rouskin_freq[mask]
		rouskin_muts = rouskin_muts[mask]
		rouskin_totals = rouskin_totals[mask]

		all_freqs += list(freq)
		all_rouskin_freqs += list(rouskin_freq)

	r = pearsonr(np.array(all_freqs), np.array(all_rouskin_freqs))[0]
	print("Pearson correlation: %f" % r)
	print("R^2 value: %f" % (r*r))

	plt.scatter(all_freqs, all_rouskin_freqs, color='black', s=2)
	plt.show()

# Get ROC curve for classifying base-paired from unpaired+accessible 
# by scanning different cutoffs for mutational frequencies
def plot_roc_auc_combined(mut_freq_file, rrna_names = ['18s', '25s']):
	bp_reacs = []
	unpaired_reacs = []

	is_base_paired_xtal_all = []
	mask_all = []
	for rrna_name in rrna_names:
		fasta_name = fasta_name_dict[rrna_name]
		is_base_paired_xtal = get_is_base_paired_xtal(rrna_name)	
		reactivities, _, _, _ = get_mut_freq(mut_freq_file, fasta_name)

		for ii, xtal_bp in enumerate(is_base_paired_xtal):
			if xtal_bp == 1:
				bp_reacs += [float(reactivities[ii])]
			if xtal_bp == 0:
				unpaired_reacs += [float(reactivities[ii])]

		mask = is_base_paired_xtal != -1
		is_base_paired_xtal = is_base_paired_xtal[mask]

		is_base_paired_xtal_all += list(is_base_paired_xtal)
		mask_all += list(mask)

	mask_all = np.array(mask_all)
	is_base_paired_xtal_all = np.array(is_base_paired_xtal_all)

	bins = np.arange(0, 1, 0.001)
	tprs = []
	fprs = []
	for cutoff in bins:
		is_base_paired_exp_all = []
		for rrna_name in rrna_names:
			fasta_name = fasta_name_dict[rrna_name]
			is_base_paired_exp = \
				get_is_base_paired_mutfrac(mut_freq_file, fasta_name, cutoff)
			is_base_paired_exp_all += list(is_base_paired_exp)
		is_base_paired_exp_all = np.array(is_base_paired_exp_all)

		is_base_paired_exp_all = is_base_paired_exp_all[mask_all]

		tp = sum(is_base_paired_exp_all * is_base_paired_xtal_all)
		tpr = tp/sum(is_base_paired_xtal_all)
		fp = sum(is_base_paired_exp_all * (1 - is_base_paired_xtal_all))
		fpr = fp/sum(1 - is_base_paired_xtal_all)
		tprs += [tpr]
		fprs += [fpr]

	plt.figure(figsize=(8, 6))
	plt.hist(bp_reacs, bins=100, alpha=0.5, label='Base-paired')
	plt.hist(unpaired_reacs, bins=100, alpha=0.5, label='Unpaired')
	plt.yscale('log')
	plt.show()

	plt.plot(fprs, tprs, color='black')
	plt.show()
	auc = metrics.auc(fprs, tprs)
	print("AUC for " + str(rrna_names) + ": " + str(auc))

	return fprs, tprs

plot_rrna_correlation()
plot_roc_auc_combined("combined_1221/rfcount/" + EXP + "_rrna_all_view.txt")
