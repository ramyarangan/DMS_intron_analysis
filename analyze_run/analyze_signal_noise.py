import numpy as np 
import sys

exp_read_lens = {'d3': 300, 'd': 300, 'nd': 300}

EXP = sys.argv[1] # n3, d3, d45, d

intron_fasta_file = "annot_files/standard_introns.fa"
mut_freq_file = "combined_1221/rfcount/" + EXP + "_all_dedup_view.txt"
stats_file = "combined_1221/run_data/" + EXP + "_stats.txt"

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
		if len(stat) > 0 and (stat[0] == 'c' or stat[0] == 'R'):
			name = stat.split('\t')[0]
			coverage = int(stat.split('\t')[1].replace('\n',''))
			coverage = coverage * exp_read_lens[experiment]/lengths_dict[name]
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

def get_acgu_mut_freq_avg_refseq(mut_freq_file, cov_dict, cov_cutoff=20):
	f = open(mut_freq_file)
	mut_freq_lines = f.readlines()
	f.close()

	totals_dict = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
	freq_dict = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

	for ii in range(int(len(mut_freq_lines)/5)):
		name = mut_freq_lines[ii*5].replace('\n', '')
		if cov_dict[name] < cov_cutoff:
			continue
		
		seq = list(mut_freq_lines[ii*5 + 1].replace('\n', ''))
		seq = np.array(seq)

		muts = mut_freq_lines[ii*5 + 2].split(',')
		muts = np.array([int(x) for x in muts])

		totals = mut_freq_lines[ii*5 + 3].split(',')
		totals = np.array([int(x) for x in totals])

		for cur_key in totals_dict.keys():
			mask = seq == cur_key
			mask = np.logical_and(mask, totals > 0)
			freq_dict[cur_key] += sum(muts[mask]/totals[mask])
			totals_dict[cur_key] += sum(mask)

	total_freq = 0
	for cur_key in totals_dict.keys():
		freq = freq_dict[cur_key]/totals_dict[cur_key]
		total_freq += freq
		print("Frequency of mutations in %s: %f" % (cur_key, freq))

	for cur_key in totals_dict.keys():
		freq = freq_dict[cur_key]/totals_dict[cur_key]
		print("Normalized mut freq in %s: %f" % (cur_key, freq/total_freq))

cov_dict = fill_coverage_dict(stats_file, EXP, intron_fasta_file)

get_acgu_mut_freq_avg_refseq(mut_freq_file, cov_dict)

