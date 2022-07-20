"""
Run statistical comparisons between intron and coding sequences for all secondary structure metrics
Make violin plots comparing distributions

Usage exameee: 
python compare_secstruct_features.py ../intron_annot/standard_introns.fa rnastructure_sherlock_1221/intron/ ../analyze_run/combined_1221/reactivity/reactivity_intron/ ../analyze_run/combined_1221/rfcount/d_all_dedup_view.txt mee_cache/intron_1221.txt ../intron_annot/standard_introns.dat ../intron_annot/coding_orfs_introns.fa rnastructure_sherlock_1221/coding_nd/ ../analyze_run/combined_1221/reactivity/reactivity_coding_nd/ ../analyze_run/combined_1221/rfcount/nd_coding_dedup_view.txt mee_cache/coding_1221.txt
"""
import sys
from matplotlib import pyplot as plt 
from scipy import stats 
import pandas as pd 
import numpy as np
import seaborn as sns 

from secstruct_util import *
from get_secstruct_features import *
from zipper_stem_feature import *
from get_gini_features import *

fasta_file = sys.argv[1] # "../intron_annot/standard_introns.fa"
secstruct_dir = sys.argv[2] # 
cov_reac_dir = sys.argv[3]
mut_freq_file = sys.argv[4] # '../analyze_run/combined_1221/rfcount/d_all_dedup_view.txt'
mee_file = sys.argv[5] # 'mee_calculations/intron_properties_600_cutoff_1221.txt'
dat_file = sys.argv[6] # "../intron_annot/standard_introns.dat"

fasta_file_coding = sys.argv[7]
secstruct_dir_coding = sys.argv[8]
cov_reac_dir_coding = sys.argv[9]
mut_freq_file_coding = sys.argv[10] # '../analyze_run/combined_1221/rfcount/nd_coding_dedup_view.txt'
mee_file_coding = sys.argv[11] # 'mee_calculations/coding_properties_600_cutoff.txt'


OVERHANG = 0
if len(sys.argv) > 12:
	OVERHANG = int(sys.argv[12])

COVERAGE_CUTOFF = 1971
snorna_introns = ['chrI:142253-142619', 'chrVII:364964-365432', 'chrXIII:499878-500151', 'chrXI:283095-283421',\
	'chrXII:856574-857057', 'chrXIII:163308-163716', 'chrXIV:721770-722302']

def get_all_bpps(names, secstruct_dir, cov_reac_dir, len_cutoff=6, do_overhang=True):
	all_bpps = []

	for ii, name in enumerate(names):
		coverage = get_cov(name, cov_reac_dir)
		if name in snorna_introns:
			continue
		if coverage > COVERAGE_CUTOFF:
			overhang = 0
			if do_overhang:
				overhang = OVERHANG
			new_bpps = get_bpps(name, secstruct_dir, len_cutoff=len_cutoff, \
				overhang_5p=overhang, overhang_3p=overhang)
			if new_bpps is not None:
				all_bpps += new_bpps

	return all_bpps

def get_all_stem_lengths(names, cov_reac_dir, secstruct_dir, do_overhang=True):
	all_stem_lens = []

	for ii, name in enumerate(names):
		coverage = get_cov(name, cov_reac_dir)
		if name in snorna_introns:
			continue
		if coverage > COVERAGE_CUTOFF:
			dotbracket = get_dotbracket(name, secstruct_dir)
			if dotbracket is None:
				continue
			bpp_matrix = get_bpp_matrix(name, secstruct_dir)
			if bpp_matrix is None:
				continue
			stemG, nt_to_node_dict = get_stem_graph(dotbracket, bpp_matrix, \
				stem_verbose=False, graph_verbose=False)
			for cur_nt in nt_to_node_dict.keys():
				cur_node = nt_to_node_dict[cur_nt]
				if cur_node.get_type() != 'stem':
					continue
				overhang = 0
				if do_overhang:
					overhang = OVERHANG
				if cur_node.get_min() < overhang or \
					cur_node.get_max() > len(dotbracket) - overhang:
					continue
				if cur_node.bpp < 0.9: 
					continue
				all_stem_lens += [len(cur_node.strand1_nts)]

	return all_stem_lens

def get_all_longest_stems(names, cov_reac_dir, secstruct_dir, do_overhang=True):
	all_longest_stems = []

	for ii, name in enumerate(names):
		coverage = get_cov(name, cov_reac_dir)
		if name in snorna_introns:
			continue
		if coverage > COVERAGE_CUTOFF:
			dotbracket = get_dotbracket(name, secstruct_dir)
			if dotbracket is None:
				continue
			bpp_matrix = get_bpp_matrix(name, secstruct_dir)
			if bpp_matrix is None:
				continue
			stemG, nt_to_node_dict = get_stem_graph(dotbracket, bpp_matrix, \
				stem_verbose=False, graph_verbose=False)
			overhang = 0
			if do_overhang:
				overhang = OVERHANG
			_, longest_stem_len = get_longest_stem(stemG, nt_to_node_dict, len(dotbracket), \
				loop_cutoff=2, bpp_cutoff=0.9, overhang_5p=overhang, \
				overhang_3p=overhang)
			if longest_stem_len == 0:
				continue

			all_longest_stems += [longest_stem_len]

	return all_longest_stems

def write_mee_to_file(fasta_file, mee_file, secstruct_dir, cov_reac_dir, \
	do_overhang=True):
	f = open(mee_file, 'w')
	
	names, _, name_seq_dict, _ = get_names_seqs_from_fasta(fasta_file)

	for ii, name in enumerate(names):
		print("Computing MEE: %d of %d\n" % (ii, len(names)))
		seq = name_seq_dict[name]
		overhang = 0
		if do_overhang:
			overhang = OVERHANG
		mee = get_mee(name, seq, secstruct_dir, cov_reac_dir, \
			overhang_5p=overhang, overhang_3p=overhang)
		f.write("%s\n" % name)
		f.write("%f\n" % mee)

	f.close()

def get_all_mee_from_file(fasta_file, mee_stats_file, min_len=-1):
	_, _, name_seq_dict, _ = get_names_seqs_from_fasta(fasta_file)
	
	f = open(mee_stats_file)
	mee_lines = f.readlines()
	f.close()

	mee_dict = {}

	for ii in range(int(len(mee_lines)/2)):
		name = mee_lines[2 * ii].replace('\n', '')
		seq_len = len(name_seq_dict[name])
		mee = float(mee_lines[2 * ii + 1])
		if mee != -1 and seq_len > min_len:
			mee_dict[name] = mee

	return mee_dict

def get_all_gini_coeffs(names, name_seq_dict, mut_freq_file, do_overhang=True):
	all_gini_coeffs = []

	for name in names:
		if name in snorna_introns:
			continue

		seq = name_seq_dict[name]
		overhang = 0
		if do_overhang:
			overhang = OVERHANG
		gini_coeffs = get_gini_windows(name, seq, mut_freq_file, \
			overhang_5p=overhang, overhang_3p=overhang)

		all_gini_coeffs += gini_coeffs

	return all_gini_coeffs

def get_all_zipper_stem_dG(names, name_seq_dict, secstruct_dir, cov_reac_dir, \
	do_second=False):
	dGs = []

	for name in names:
		coverage = get_cov(name, cov_reac_dir)

		if coverage < COVERAGE_CUTOFF:
			continue

		mfe = get_dotbracket(name, secstruct_dir)
		bpp_matrix = get_bpp_matrix(name, secstruct_dir)
		bp = get_bp_loc(name, dat_file)

		has_zipper_stem = False
		if (mfe is not None) and (bpp_matrix is not None) and (bp > 0):
			if do_second:
				has_zipper_stem, stem, best_dG = get_best_zipper_stem(bp, name_seq_dict[name], \
					mfe, bpp_matrix, min_num_bp=6, do_second=True, \
					overhang_5p=OVERHANG, overhang_3p=OVERHANG)
			else:
				has_zipper_stem, stem, best_dG = get_best_zipper_stem(bp, name_seq_dict[name], \
					mfe, bpp_matrix, min_num_bp=6, overhang_5p=OVERHANG, \
					overhang_3p=OVERHANG)

		dG = ""
		if has_zipper_stem:
			dG = str(best_dG)

		dGs += [dG]

	return dGs

def plot_violin_plot(intron_vals, coding_vals, plt_title, savefig=False, filename='', do_inner=True):
	df = pd.DataFrame(columns=['Intron', 'Coding'])

	intron_added_vals = [np.nan] * max(len(intron_vals), len(coding_vals))
	intron_added_vals[:len(intron_vals)] = intron_vals 
	df = df.assign(Intron = intron_added_vals)
	coding_added_vals = [np.nan] * max(len(intron_vals), len(coding_vals))
	coding_added_vals[:len(coding_vals)] = coding_vals
	df = df.assign(Coding = coding_added_vals)
	plt.figure(figsize=(4, 3))
	my_pal = {"Intron": "purple", "Coding": "darkorange"}
	ax = sns.violinplot(data=df, palette=my_pal, cut=0)
	if do_inner:
		ax = sns.violinplot(data=df, inner="points", palette=my_pal, cut=0)
	plt.title(plt_title)
	if savefig:
		if filename == '':
			raise RuntimeError("Must provide figure filename")
		plt.savefig(filename + '.png', format='png', dpi=300)
	else:
		plt.show()

def compare_intron_coding_vals(intron_vals, coding_vals, plt_title, \
	do_plot=True, savefig=False, filename=""):
	stat, pval = stats.mannwhitneyu(intron_vals, coding_vals, alternative='greater')
	print(pval)

	print("Number of intron stems passing cutoff: %d" % len(intron_vals))
	print("Number of coding stems passing cutoff: %d" % len(coding_vals))
	plot_violin_plot(intron_vals, coding_vals, plt_title, do_inner=False)
	if savefig:
		plot_violin_plot(intron_vals, coding_vals, plt_title, do_inner=False, 
			savefig=savefig, filename=filename)

def longest_stem_compare():
	names, _, _, _ = get_names_seqs_from_fasta(fasta_file)
	intron_stem_lens = get_all_longest_stems(names, cov_reac_dir, secstruct_dir)

	names, _, _, _ = get_names_seqs_from_fasta(fasta_file_coding)
	coding_stem_lens = get_all_longest_stems(names, cov_reac_dir_coding, \
		secstruct_dir_coding, do_overhang=False)

	print("Mann Whitney U rank test p value comparing intron / coding longest stem lengths:")
	compare_intron_coding_vals(intron_stem_lens, coding_stem_lens, "Length of longest stem", 
		filename='../figures/longest_stem_len_compare', savefig=True)

def stem_len_compare():
	names, _, _, _ = get_names_seqs_from_fasta(fasta_file)
	intron_stem_lens = get_all_stem_lengths(names, cov_reac_dir, secstruct_dir)

	names, _, _, _ = get_names_seqs_from_fasta(fasta_file_coding)
	coding_stem_lens = get_all_stem_lengths(names, cov_reac_dir_coding, \
		secstruct_dir_coding, do_overhang=False)

	print("Mann Whitney U rank test p value comparing intron / coding stem lengths:")
	compare_intron_coding_vals(intron_stem_lens, coding_stem_lens, "Stem lengths",
		filename='../figures/stem_len_compare', savefig=True)

def bpp_compare(len_cutoff=6):
	names, _, _, _ = get_names_seqs_from_fasta(fasta_file)
	intron_bpps = get_all_bpps(names, secstruct_dir, cov_reac_dir, len_cutoff=len_cutoff)

	names, _, _, _ = get_names_seqs_from_fasta(fasta_file_coding)
	coding_bpps = get_all_bpps(names, secstruct_dir_coding, cov_reac_dir_coding, \
		len_cutoff=len_cutoff, do_overhang=False)

	print("Mann Whitney U rank test p value comparing intron / coding base-pair probability in stems:")
	compare_intron_coding_vals(intron_bpps, coding_bpps, "BPP comparison", 
		filename='../figures/bpp_compare', savefig=True)

def gini_compare():
	names, _, name_seq_dict, _ = get_names_seqs_from_fasta(fasta_file)
	intron_gini = get_all_gini_coeffs(names, name_seq_dict, mut_freq_file)

	names, _, name_seq_dict, _ = get_names_seqs_from_fasta(fasta_file_coding)
	coding_gini = get_all_gini_coeffs(names, name_seq_dict, mut_freq_file_coding, \
		do_overhang=False)

	print("Mann Whitney U rank test p value comparing intron / coding Gini coeffs:")
	compare_intron_coding_vals(intron_gini, coding_gini, "Gini Coefficient comparison", 
		filename='../figures/Gini', savefig=True)

def mee_compare():
	if not os.path.exists(mee_file):
		write_mee_to_file(fasta_file, mee_file, secstruct_dir, cov_reac_dir)
	if not os.path.exists(mee_file_coding):
		write_mee_to_file(fasta_file_coding, mee_file_coding, \
			secstruct_dir_coding, cov_reac_dir_coding, do_overhang=False)

	intron_mee_dict = get_all_mee_from_file(fasta_file, mee_file, min_len=0)
	coding_mee_dict = get_all_mee_from_file(fasta_file_coding, mee_file_coding, min_len=0)
	mee_normalized_intron = list(intron_mee_dict.values())
	mee_normalized_coding = list(coding_mee_dict.values())

	print("Mann Whitney U rank test p value comparing normalized MEE:")
	compare_intron_coding_vals(mee_normalized_intron, mee_normalized_coding, \
		"Normalized Maximum Extrusion from Ends", filename='../figures/MEE_normalized_violinplot', savefig=True)

def zipper_stem_stats():
	names, _, name_seq_dict, _ = get_names_seqs_from_fasta(fasta_file)

	dGs = get_all_zipper_stem_dG(names, name_seq_dict, secstruct_dir, cov_reac_dir)
	dGs = np.array(dGs)
	num_zipper_stems = np.sum(dGs != "")

	print("Number of zipper stems: %d\n" % num_zipper_stems)

	dGs = get_all_zipper_stem_dG(names, name_seq_dict, secstruct_dir, cov_reac_dir, do_second=True)
	dGs = np.array(dGs)
	num_end_stems = np.sum(dGs != "")

	print("Number of end stems: %d\n" % num_end_stems)

longest_stem_compare()
# stem_len_compare()
bpp_compare()
gini_compare()
zipper_stem_stats()
mee_compare()
