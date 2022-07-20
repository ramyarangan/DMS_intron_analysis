"""
Compare native secondary structures to DMS-guided secondary structure predictions
for control RNA's, computing confusion matrix (TP, FP, FN) for stem predictions.
"""
import os
import numpy as np
from matplotlib import pyplot as plt
from arnie import mfe 

BASE_PAIR_SYMBOLS = {"}": "{", ")": "(", "]": "[", ">": "<"}

# Get the sequence tags and sequences from a fasta file in 
# a couple data structures (lists and a dictionary)
def get_names_seqs_from_fasta(fasta_file):
	f = open(fasta_file)
	lines = f.readlines()
	f.close()

	names = []
	seqs = []
	name_seq_dict = {}
	cur_seq = ''
	for line in lines:
		if line[0] == '>':
			if cur_seq != '':
				seqs += [cur_seq]
				name_seq_dict[names[-1]] = cur_seq
			cur_seq = ''
			name = line.split(' ')[0].replace('>', '').replace('\n', '')
			names += [name]
		else:
			cur_seq += line.replace('\n', '')
	seqs += [cur_seq]
	name_seq_dict[names[-1]] = cur_seq
	return names, seqs, name_seq_dict

# Get secondary structure from file
def get_dotbracket(name, secstruct_dir):
	dotbracket_filename = secstruct_dir + name + '_secstruct.txt'
	if not os.path.exists(dotbracket_filename):
		return None

	f = open(dotbracket_filename)
	lines = f.readlines()
	f.close()
	return lines[0].replace('\n', '')

# Get base-pair probability matrix from file
def get_bpp_matrix(name, secstruct_dir):
	bpp_matrix_filename = secstruct_dir + name + '_bpp.csv'
	if not os.path.exists(bpp_matrix_filename):
		return None 

	f = open(bpp_matrix_filename)
	lines = f.readlines()
	f.close()
	
	bpp_matrix = []
	for line in lines: 
		bpp_matrix_line = line.split(',')
		bpp_matrix_line = [float(x) for x in bpp_matrix_line]
		bpp_matrix += [bpp_matrix_line]
	
	return bpp_matrix

# Get minimum free energy structure prediction
# This will be compared to the native structure along with
# DMS-guided prediction
def get_mfe(name, fasta_file, package="vienna"):
	_, _, name_to_seq_dict = get_names_seqs_from_fasta(fasta_file)
	seq = name_to_seq_dict[name]
	return mfe.mfe(seq, package=package)

# Stem class: 5' and 3' strand nucleotides along with basic functions
class Stem:
	def __init__(self, strand1_nts, strand2_nts):
		if len(strand1_nts) != len(strand2_nts):
			raise RuntimeError("Stem needs to have equal # nts in strands")

		self.strand1_nts = strand1_nts
		self.strand2_nts = strand2_nts

	def get_bpp(self, bpp_matrix):
		bpp = 0
		for ii, nt1 in enumerate(self.strand1_nts):
			nt2 = self.strand2_nts[ii]
			bpp = max(bpp, bpp_matrix[nt1][nt2])
		return bpp

	def len(self):
		return len(self.strand1_nts)

	def is_in(self, start_idx, end_idx):
		min_within = (max(self.strand1_nts) < end_idx) and \
			(max(self.strand1_nts) > start_idx)
		max_within = (min(self.strand2_nts) < end_idx) and \
			(min(self.strand2_nts) > start_idx)
		return min_within or max_within

	def contains(self, nt):
		return (nt in self.strand1_nts) or (nt in self.strand2_nts)

	def __str__(self):
		print_str = "Stem of length %d containing base-pairs:\n" % self.len()
		for ii, n1 in enumerate(self.strand1_nts):
			print_str += "%d %d\n" % (n1, self.strand2_nts[ii])
		return print_str

# Get all stems from the secondary structure in dot-bracket notation
# Allows for bulges of size at most max_bulge_cnt within a stem
# Returns a list of stem objects along with a dictionary from nucleotide
# index to the stem containing it
def get_stems(dotbracket, max_bulge_cnt=0):
	base1_stacks = {'{': [], '(': [], '[': [], '<': []}
	cur_stems = {'{': [], '(': [], '[': [], '<': []}
	
	stem_list = []

	end_bulge_cnt = 0
	for ii, curchar in enumerate(dotbracket):
		in_stem = False

		if curchar in BASE_PAIR_SYMBOLS.keys():
			end_bulge_cnt = 0
			cur_start = BASE_PAIR_SYMBOLS[curchar]

			cur_stem = cur_stems[cur_start]
			if len(cur_stem) > 0 and \
				((base1_stacks[cur_start][-1] <= cur_stem[-1][0] - 1) and 
				(base1_stacks[cur_start][-1] >= cur_stem[-1][0] - 1 - max_bulge_cnt)):
				in_stem = True

		if curchar == "." or curchar == "-":
			if end_bulge_cnt < max_bulge_cnt:
				in_stem = True
			end_bulge_cnt += 1

		if not in_stem:
			for cur_start, cur_stem in cur_stems.items():
				if len(cur_stem) > 0:
					strand1_nts = [x[0] for x in cur_stem]
					strand2_nts = [x[1] for x in cur_stem]
					stem_list += [Stem(strand1_nts, strand2_nts)]
					cur_stems[cur_start] = []
					end_bulge_cnt = 0

		if curchar in base1_stacks.keys():
			base1_stacks[curchar] += [ii]

		if curchar in BASE_PAIR_SYMBOLS.keys():
			cur_start = BASE_PAIR_SYMBOLS[curchar]
			cur_stems[cur_start] += [(base1_stacks[cur_start][-1], ii)]
			base1_stacks[cur_start] = base1_stacks[cur_start][:-1]

	for cur_start, cur_stem in cur_stems.items():
		if len(cur_stem) > 0:
			strand1_nts = [x[0] for x in cur_stem]
			strand2_nts = [x[1] for x in cur_stem]
			stem_list += [Stem(strand1_nts, strand2_nts)]

	nt_to_stem_dict = {}

	for stem in stem_list:
		for nt in stem.strand1_nts:
			nt_to_stem_dict[nt] = stem
		for nt in stem.strand2_nts:
			nt_to_stem_dict[nt] = stem

	return stem_list, nt_to_stem_dict

# Get all base pairs in a secondary structure, returns
# as a dictionary from a nucleotide to its base-pairing partner 
def get_base_pairs(dotbracket):
	basepair_dict = {}
	base1_stacks = {'{': [], '(': [], '[': [], '<': []}

	for ii, curchar in enumerate(dotbracket):
		if curchar in base1_stacks.keys():
			base1_stacks[curchar] += [ii]

		if curchar in BASE_PAIR_SYMBOLS.keys():
			start_symb = BASE_PAIR_SYMBOLS[curchar]
			bp_partner = base1_stacks[start_symb][-1]
			base1_stacks[start_symb] = base1_stacks[start_symb][:-1]
			basepair_dict[ii] = bp_partner
			basepair_dict[bp_partner] = ii

	return basepair_dict

# Get the percent of base-pairs in the stem that are matched 
# in the base-pair dictionary
def get_perc_match(stem, bp_dict):
	strand1_nts = stem.strand1_nts
	strand2_nts = stem.strand2_nts

	num_matched = 0
	matched_nt = -1
	for ii, nt in enumerate(strand1_nts):
		if nt in bp_dict.keys() and bp_dict[nt] == strand2_nts[ii]:
			num_matched += 1
			matched_nt = nt

	return num_matched/stem.len(), matched_nt

# Compares two secondary structures and identifies true positive, 
# false positive, and false negative rates
def get_conf_matrix_from_stem_list(dms_struct, native_struct, bpp_matrix, \
	limits, min_stem_len, perc_match_cutoff, bootstrap_cutoff):
	dms_stem_list, _ = get_stems(dms_struct)
	bp_dict = get_base_pairs(native_struct)
	native_stem_list, _ = get_stems(native_struct)

	tp = 0
	fp = 0
	# Positions marking stems in the native structure that have found a match 
	# with the DMS MaP-seq structure
	matched_nts = []
	fn_nts = []
	for stem in dms_stem_list:
		bpp = stem.get_bpp(bpp_matrix)

		if stem.len() > min_stem_len:
			# Only evaluate stems that are within the limits
			# of the control structure
			stem_passes = True
			if limits[0] != -1:
				stem_passes = stem.is_in(limits[0], limits[1])
			if not stem_passes:
				continue

			perc_match, nt = get_perc_match(stem, bp_dict)

			if bpp < bootstrap_cutoff:
				# A stem that matches the native structure but has
				# low confidence estimate is a false negative
				if perc_match >= perc_match_cutoff:
					# Separate FN into classes using this
					fn_nts += [nt]
				continue

			# Only predicted stems with confidence estimates above the
			# bootstrap cutoff are positive predictions
			if perc_match >= perc_match_cutoff:
				matched_nts += [nt]
				tp += 1
			else:
				fp += 1

	fn = 0
	fn_low_boostrap = 0

	# Go through stems in the native structure to find false negatives
	for stem in native_stem_list:
		stem_in_dms = False
		stem_in_struct_low_bootstrap = False
		for nt in matched_nts:
			if stem.contains(nt):
				stem_in_dms = True
				break
		for nt in fn_nts:
			# Also tabulate false negatives that were in a predicted stem, 
			# but the stem had low helix confidence estimate
			if stem.contains(nt):
				stem_in_struct_low_bootstrap = True
				break

		if (stem.len() > min_stem_len) and (not stem_in_dms):
			fn += 1
	
		if (stem.len() > min_stem_len) and stem_in_struct_low_bootstrap:
			fn_low_boostrap += 1

	return tp, fp, fn, 0, fn_low_boostrap


# We require min_stem_len base-pairs in a helix of 70% bootstrapping probability
def get_confusion_matrix(name, dms_secstruct_dir, native_secstruct_dir, limits, \
	min_stem_len=4, bootstrap_cutoff=0.7, perc_match_cutoff=0.5, verbose=False, \
	fasta_file="", do_mfe=False):
	dms_struct = get_dotbracket(name, dms_secstruct_dir)
	if do_mfe:
		dms_struct = get_mfe(name, fasta_file)
	native_struct = get_dotbracket(name, native_secstruct_dir)
	if verbose:
		print(len(dms_struct))

	bpp_matrix = get_bpp_matrix(name, dms_secstruct_dir)
	bpp_matrix = np.array(bpp_matrix)

	tp, fp, fn, tn, fn_low_boostrap = get_conf_matrix_from_stem_list(dms_struct, \
		native_struct, bpp_matrix, limits, min_stem_len, perc_match_cutoff, bootstrap_cutoff)

	return tp, fp, fn, tn, fn_low_boostrap

# Get confusion matrix for each control RNA and make a combined matrix
def get_total_confusion_matrix(fasta_file="", do_mfe=False, bootstrap_cutoff=0.7, \
	print_summary=False, verbose=False):
	control_names = ["RDN18-1", "RDN5-1", "RDN58-1", "SNR7-L", "SNR19", "HAC1", \
		"ASH1", "RPS28B", "SFT2", "TRR4", "TRT2", "IMT4"]
	limits = {}
	for control_name in control_names:
		limits[control_name] = (-1, -1)
	limits["HAC1"] = (150, 232)
	limits["ASH1"] = (150, 206)
	limits["RPS28B"] = (150, 201)
	limits["SFT2"] = (150, 182)
	mrna_controls = ["HAC1", "ASH1", "RPS28B", "SFT2"]

	total_conf_matrix = np.array([0] * 5)
	total_conf_matrix_mrna = np.array([0] * 5)
	for control_name in control_names:
		confusion_matrix = get_confusion_matrix(control_name, "control_secstruct_bpps/", 
			"native_secstructs/", limits[control_name], fasta_file=fasta_file, \
			do_mfe=do_mfe, bootstrap_cutoff=bootstrap_cutoff, verbose=verbose)
		if verbose:
			print(control_name)
			print(confusion_matrix)
		total_conf_matrix += confusion_matrix
		if control_name in mrna_controls:
			total_conf_matrix_mrna += confusion_matrix

	if print_summary:
		print("Total confusion matrix")
		print(total_conf_matrix)
	if verbose:
		print("Total confusion matrix mRNA")
		print(total_conf_matrix_mrna)

	return total_conf_matrix

def get_accuracy(confusion_matrix):
	(tp, fp, fn, tn) = tuple(confusion_matrix)
	return (tn + tp)/(tn + fp + tp + fn)

def get_precision(confusion_matrix):
	(tp, fp, fn, tn) = tuple(confusion_matrix)
	return tp/(tp + fp)

def get_recall(confusion_matrix):
	(tp, fp, fn, tn) = tuple(confusion_matrix)
	return tp/(tp + fn)

def get_f1_score(confusion_matrix):
	precision = get_precision(confusion_matrix)
	recall = get_recall(confusion_matrix)
	return 2 * (precision * recall)/(precision + recall)

# Plot accuracy, precision, and recall over a range of
# helix confidence estimate cutoff scores, adding in the values for 
# Vienna MFE predictions as red lines
def plot_metrics(bootstrap_range):
	accuracies = []
	precisions = []
	recalls = []
	f1_scores = []

	for cutoff in bootstrap_range:
		conf_matrix = get_total_confusion_matrix(bootstrap_cutoff=cutoff)
		accuracies += [get_accuracy(conf_matrix[0:4])]
		precisions += [get_precision(conf_matrix[0:4])]
		recalls += [get_recall(conf_matrix[0:4])]
		f1_scores += [get_f1_score(conf_matrix[0:4])]

	fig, axs = plt.subplots(1, 3)
	fig.set_size_inches(18, 4)

	axs[0].plot(bootstrap_range, precisions, color='black', linewidth=1.2)
	axs[0].axvline(0.7, linestyle='--', color='black', linewidth=0.6)
	axs[0].axhline(0.533, linestyle='--', color='red', linewidth=0.6)
	axs[0].set_title("PPV for stem prediction")
	axs[0].set_xlabel("Helix confidence estimate cutoff")
	axs[0].set_ylabel("PPV")
	axs[0].set_ylim((0.5, 0.9))

	axs[1].plot(bootstrap_range, recalls, color='black', linewidth=1.2)
	axs[1].axvline(0.7, linestyle='--', color='black', linewidth=0.6)
	axs[1].axhline(0.75, linestyle='--', color='red', linewidth=0.6)
	axs[1].set_title("Sensitivity for stem prediction")
	axs[1].set_xlabel("Helix confidence estimate cutoff")
	axs[1].set_ylabel("Sensitivity")
	axs[1].set_ylim((0.5, 0.9))

	# axs[1, 0].plot(bootstrap_range, accuracies, color='black')
	# axs[1, 0].axvline(0.7, linestyle='--', color='red')
	# axs[1, 0].set_title("Accuracy")
	# axs[1, 0].set_xlabel("Helix confidence estimate cutoff")
	# axs[1, 0].set_ylabel("Accuracy for stem prediction")

	axs[2].plot(bootstrap_range, f1_scores, color='black', linewidth=1.2)
	axs[2].axvline(0.7, linestyle='--', color='black', linewidth=0.6)
	axs[2].axhline(0.62, linestyle='--', color='red', linewidth=0.6)
	axs[2].set_title("F1 score for stem prediction")
	axs[2].set_xlabel("Helix confidence estimate cutoff")
	axs[2].set_ylabel("F1 score")
	axs[2].set_ylim((0.5, 0.9))

	print(precisions)
	print(bootstrap_range)
	print(f1_scores)
	# plt.show()
	plt.savefig("../figures/control_ppv_sens_f1.png", format='png', dpi=300)

# Get precision, recall, F1 score for Vienna MFE predictions
print("Vienna MFE prediction")
conf_matrix = get_total_confusion_matrix(fasta_file="../intron_annot/control_RNAs_predicted.fa", do_mfe=True, \
	print_summary=True, bootstrap_cutoff=0)
print("Precision: %f" % get_precision(conf_matrix[0:4]))
print("Recall: %f" % get_recall(conf_matrix[0:4]))
print("F1 Score: %f" % get_f1_score(conf_matrix[0:4]))

# Get confusion matrix for raw DMS-guided structure prediction
print("DMS-guided RNAstructure, no helix confidence estimate cutoff")
conf_matrix = get_total_confusion_matrix(print_summary=True, bootstrap_cutoff=0)
print("Precision: %f" % get_precision(conf_matrix[0:4]))
print("Recall: %f" % get_recall(conf_matrix[0:4]))
print("F1 Score: %f" % get_f1_score(conf_matrix[0:4]))

# Get confusion matrix for DMS-guided structure prediction with helix confidence estimate cutoff
print("DMS-guided RNAstructure, helix confidence estimate cutoff 0.7")
conf_matrix = get_total_confusion_matrix(print_summary=True)
print("Precision: %f" % get_precision(conf_matrix[0:4]))
print("Recall: %f" % get_recall(conf_matrix[0:4]))
print("F1 score: %f" % get_f1_score(conf_matrix[0:4]))

# Plot precision, recall, F1 score over a range of helix confidence estimate cutoffs
bootstrap_range = np.arange(0, 1, 0.1)
plot_metrics(bootstrap_range)

