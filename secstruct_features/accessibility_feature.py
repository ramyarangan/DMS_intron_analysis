"""
Compute average accessibility at 5' splice site, branchpoint, and 3' splice site
"""
from secstruct_util import *
import numpy as np
from scipy.stats import pearsonr
from matplotlib import pyplot as plt 

# 1971 is the requirement per replicate, so we ask for twice this total
# Note - coverage in the extended alignment file is from both replicates
COV_CUTOFF = 1971 * 2

# Parameters for calculating accessibility of key sequences
FIVESS_LEN = 6
BP_OFFSET = -6
BP_LEN = 8
THREESS_OFFSET = -4
THREESS_LEN = 4

# Gets the base data tagline from an extended window
# that includes this intron tag
def get_extended_window_tag_ss_loc(tag, extended_base_data):
	chr_num = tag.split(":")[0]
	chr_start = tag.split(":")[1].split("-")[0]
	chr_end = tag.split(":")[1].split("-")[1]

	f = open(extended_base_data)
	basedata_lines = f.readlines()
	f.close()

	for ii in range(int(len(basedata_lines)/2)):
		tagline = basedata_lines[ii * 2]

		basedata_items = tagline.split(' ')[1].split('\t')
		cur_chr_num = basedata_items[0]
		cur_chr_start = basedata_items[1]
		cur_chr_end = basedata_items[2]

		if chr_num != cur_chr_num:
			continue

		if (int(cur_chr_start) <= int(chr_start)) and \
			(int(cur_chr_end) >= int(chr_end)):
			bp_pos = int(tagline.split(' ')[0])
			start_offset = int(basedata_items[6])
			end_offset = int(basedata_items[7])
			cur_tag = cur_chr_num + ":" + cur_chr_start + "-" + cur_chr_end + \
				"(" + basedata_items[5] + ")"
			return cur_tag, bp_pos, start_offset, end_offset

	return "", -1, -1, -1

def get_reac_bp_from_tag(tag, extended_base_data, reac_dir):
	ext_tag, bp_pos, start_offset, end_offset = \
		get_extended_window_tag_ss_loc(tag, extended_base_data)
	
	reac, cov = get_reactivities(ext_tag, reac_dir)
	
	if cov < COV_CUTOFF:
		return -1, -1, -1

	reac = np.array(reac)

	fivess_reac = reac[start_offset:(start_offset + FIVESS_LEN)]
	bp_reac = reac[(bp_pos+BP_OFFSET):(bp_pos+BP_OFFSET+BP_LEN)]
	end_pos = len(reac)-end_offset-1
	threess_reac = reac[(end_pos+THREESS_OFFSET):(end_pos+THREESS_OFFSET+THREESS_LEN)]

	has_nonnan = True
	if sum(fivess_reac[~np.isnan(fivess_reac)]) == 0:
		has_nonnan = False
	if sum(bp_reac[~np.isnan(bp_reac)]) == 0:
		has_nonnan = False
	if sum(threess_reac[~np.isnan(threess_reac)]) == 0:
		has_nonnan = False

	if not has_nonnan:
		return -1, -1, -1

	fivess_avg_reac = np.nanmean(fivess_reac)
	bp_avg_reac = np.nanmean(bp_reac)
	threess_avg_reac = np.nanmean(threess_reac)

	return fivess_avg_reac, bp_avg_reac, threess_avg_reac

def plot_correlation(vals1, vals2, plt_title):
	r = pearsonr(np.array(vals1), np.array(vals2))
	print(r)

	plt.scatter(vals1, vals2, color='black', s=0.8)
	plt.title(plt_title)
	plt.show()

def get_correlations():
	fasta_file = "../intron_annot/standard_introns.fa"
	extended_dat_file = "../intron_annot/standard_allsize_extend50_baseinfo.dat"
	extended_reac_dir = "../analyze_run/combined_1221/reactivity/reactivity_allsize_extend50/"
	ri_fractions_file = "../analyze_run/ri_fractions.txt"

	names, _, _, _ = get_names_seqs_from_fasta(fasta_file)
	pladb_ri_dict, no_pladb_ri_dict = get_ri_fractions(ri_fractions_file)

	fivess_reacs = []
	bp_reacs = []
	threess_reacs = []
	pladb_ri = []
	pladb_ri_ratio = []

	for name in names:
		fivess_reac, bp_reac, threess_reac = \
			get_reac_bp_from_tag(name, extended_dat_file, extended_reac_dir)
		if fivess_reac == -1:
			continue

		if name not in pladb_ri_dict.keys():
			continue

		fivess_reacs += [fivess_reac]
		bp_reacs += [bp_reac]
		threess_reacs += [threess_reac]

		pladb_ri += [pladb_ri_dict[name]]
		pladb_ri_ratio += [pladb_ri_dict[name]/no_pladb_ri_dict[name]]

	plot_correlation(fivess_reacs, bp_reacs, "5'SS vs BP accessibility")
	plot_correlation(fivess_reacs, threess_reacs, "5'SS vs 3'SS accessibility")
	plot_correlation(fivess_reacs, pladb_ri, "5'SS vs pladB RI fraction")
	plot_correlation(fivess_reacs, pladb_ri_ratio, "5'SS vs pladB RI ratio")

# get_correlations()


