import numpy as np 
import sys
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from matplotlib import pyplot as plt 

exp_read_counts = {'d3': 557428287, 'd45': 380745898, 'd': 938174185, 'n3': 9871773, 'nd': 169742497, 'rouskin': 296022844}
exp_read_lens = {'d3': 300, 'd45': 300, 'd': 300, 'nd': 300, 'n3': 300, 'rouskin': 50, 'd3_0.25': 300, 'd45_0.25': 300}

EXP1 = sys.argv[1] # d3
EXP2 = sys.argv[2] # d45

intron_fasta_file = "annot_files/standard_introns.fa"

mut_freq_file_1 = "combined/rfcount/" + EXP1 + "_all_dedup_view.txt"
stats_file_1 = "combined/run_data/" + EXP1 + "_stats.txt"

mut_freq_file_2 = "combined/rfcount/" + EXP2 + "_all_dedup_view.txt"
stats_file_2 = "combined/run_data/" + EXP2 + "_stats.txt"

# Get dictionary of name of construct to length of each construct from a fasta file
def fill_length_dict(filename):
	f = open(filename)
	seqs = f.readlines()
	f.close()

	len_dict = {}
	tag = ""
	cur_seq = ""
	for seq in seqs:
		if len(seq) > 0 and seq[0] == ">":
			if tag != "":
				len_dict[tag] = len(cur_seq)
			tag = seq.split(" ")[0].replace(">", "").replace('\n', '')
			cur_seq = ""
		else:
			cur_seq += seq.replace("\n", "")
	if tag != "":
		len_dict[tag] = len(cur_seq)	
	return len_dict

# RPKM mode: RPKM normalized using the construct lengths from the fasta file and experiment read counts
# Not RPKM mode: normalize to construct lengths from the fasta file, multiply by read length. This 
#               is the expected coverage per position in the construct.
def get_norm_cov(stats, lengths_dict, experiment, do_rpkm=True, do_rpk=False):
	cov_dict = {}
	for stat in stats:
		if len(stat) > 0 and stat[0] == 'c':
			name = stat.split('\t')[0]
			coverage = int(stat.split('\t')[1].replace('\n',''))
			coverage = coverage * 1000/lengths_dict[name]
			if do_rpkm:
				permil_factor = exp_read_counts[experiment]/1000000
				coverage = coverage/permil_factor
			if not do_rpkm and not do_rpk:
				coverage = coverage * exp_read_lens[experiment]
			cov_dict[name] = coverage
	return cov_dict

# Fill coverage dictionary with values from the RNAframework stats files
# RNAframework stats files is the number of reads that map onto each sequence
def fill_coverage_dict(intron_filename, experiment, fasta_file, do_rpkm=False, do_rpk=False):
	f = open(intron_filename)
	intron_stats = f.readlines()
	f.close()

	intron_lens = fill_length_dict(fasta_file)

	intron_cov_dict = get_norm_cov(intron_stats, intron_lens, experiment, do_rpkm, do_rpk)

	return intron_cov_dict

def plot_reactivity_correlation(mut_freq_file_1, mut_freq_file_2, \
	cov_dict_1, cov_dict_2, cov_cutoff=10):
	f = open(mut_freq_file_1)
	mut_freq_lines_1 = f.readlines()
	f.close()

	f = open(mut_freq_file_2)
	mut_freq_lines_2 = f.readlines()
	f.close()

	mut_freqs_1 = []
	mut_freqs_2 = []
	for ii in range(int(len(mut_freq_lines_1)/5)):
		name1 = mut_freq_lines_1[ii*5].replace('\n', '')
		name2 = mut_freq_lines_2[ii*5].replace('\n', '')
		if name1 != name2:
			raise RuntimeError

		if (cov_dict_1[name1] + cov_dict_2[name2])/2 < cov_cutoff:
			continue
		
		seq = list(mut_freq_lines_1[ii*5 + 1].replace('\n', ''))
		seq = np.array(seq)

		muts = mut_freq_lines_1[ii*5 + 2].split(',')
		muts = np.array([int(x) for x in muts])

		totals = mut_freq_lines_1[ii*5 + 3].split(',')
		totals = np.array([int(x) for x in totals])

		freqs1 = muts/totals

		muts = mut_freq_lines_2[ii*5 + 2].split(',')
		muts = np.array([int(x) for x in muts])

		totals = mut_freq_lines_2[ii*5 + 3].split(',')
		totals = np.array([int(x) for x in totals])

		freqs2 = muts/totals

		for mod_base in ['A', 'C']:
			mask = seq == mod_base
			mut_freqs_1 += list(freqs1[mask])
			mut_freqs_2 += list(freqs2[mask])
	
	r = pearsonr(np.array(mut_freqs_1), np.array(mut_freqs_2))[0]
	print(r)
	print(r*r)

	plt.scatter(mut_freqs_1, mut_freqs_2)
	plt.show()

def make_reactivity_correlation_graph(mut_freq_file_1, mut_freq_file_2, \
	cov_dict_1, cov_dict_2):
	f = open(mut_freq_file_1)
	mut_freq_lines_1 = f.readlines()
	f.close()

	f = open(mut_freq_file_2)
	mut_freq_lines_2 = f.readlines()
	f.close()

	cutoffs_list = np.arange(5, 200, 5)
	num_passing_list = []
	r2_list = []
	r2_20 = 0
	r2_40 = 0
	num_passing_20 = 0
	num_passing_40 = 0
	for cutoff in cutoffs_list:
		mut_freqs_1 = []
		mut_freqs_2 = []
		num_passing = 0
		for ii in range(int(len(mut_freq_lines_1)/5)):
			name1 = mut_freq_lines_1[ii*5].replace('\n', '')
			name2 = mut_freq_lines_2[ii*5].replace('\n', '')
			if name1 != name2:
				raise RuntimeError

			if (cov_dict_1[name1] + cov_dict_2[name2])/2 < cutoff:
				continue
			num_passing += 1
			seq = list(mut_freq_lines_1[ii*5 + 1].replace('\n', ''))
			seq = np.array(seq)

			muts = mut_freq_lines_1[ii*5 + 2].split(',')
			muts = np.array([int(x) for x in muts])

			totals = mut_freq_lines_1[ii*5 + 3].split(',')
			totals = np.array([int(x) for x in totals])

			freqs1 = muts/totals

			muts = mut_freq_lines_2[ii*5 + 2].split(',')
			muts = np.array([int(x) for x in muts])

			totals = mut_freq_lines_2[ii*5 + 3].split(',')
			totals = np.array([int(x) for x in totals])

			freqs2 = muts/totals

			for mod_base in ['A', 'C']:
				mask = seq == mod_base
				mut_freqs_1 += list(freqs1[mask])
				mut_freqs_2 += list(freqs2[mask])
		
		r = pearsonr(np.array(mut_freqs_1), np.array(mut_freqs_2))[0]
		if cutoff == 20:
			r2_20 = r
			num_passing_20 = num_passing
		if cutoff == 40:
			r2_40 = r
			num_passing_40 = num_passing
		r2_list += [r]
		num_passing_list += [num_passing]

	plt.plot(cutoffs_list, r2_list, color='black')
	plt.axhline(y=r2_20, color='r', linestyle='--')
	plt.axhline(y=r2_40, color='r', linestyle='--')
	plt.show()

	plt.plot(cutoffs_list, num_passing_list, color='black')
	plt.axhline(y=num_passing_20, color='r', linestyle='--')
	plt.axhline(y=num_passing_40, color='r', linestyle='--')
	plt.show()

cov_dict_1 = fill_coverage_dict(stats_file_1, EXP1, intron_fasta_file, do_rpkm=True)
cov_dict_2 = fill_coverage_dict(stats_file_2, EXP2, intron_fasta_file, do_rpkm=True)
# plot_reactivity_correlation(mut_freq_file_1, mut_freq_file_2, cov_dict_1, cov_dict_2)
make_reactivity_correlation_graph(mut_freq_file_1, mut_freq_file_2, cov_dict_1, cov_dict_2)
