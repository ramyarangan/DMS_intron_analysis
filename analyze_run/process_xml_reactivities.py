"""
Given xml files from RNAFramework's containing reactivity values, write a .txt file with 
reactivity values and read coverage across the sequence for every sequence in the fasta file. 

Last modified: 1/12/22
Example usage: 
python process_xml_reactivities.py ../intron_annot/standard_allsize_extend50.fa combined_1221/rfnorm_reactivity/rfnorm_d_allsize_extend50/ combined_1221/run_data/d_allsize_extend50_stats.txt combined_1221/reactivity/reactivity_allsize_extend50
"""
import os
import sys
import numpy as np
import xml.etree.ElementTree as ET

READ_LEN = 300

fasta_file = sys.argv[1]
reac_dir = sys.argv[2]
stats_file = sys.argv[3]
outdir = sys.argv[4]

window_start = 0
window_end = 0
if len(sys.argv) > 5:
	window_start = int(sys.argv[5])
	window_end = int(sys.argv[6])

# Gets sequences in FASTQ file into a couple different data structures
def get_names_seqs(fasta_file):
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

def get_norm_cov(stats, lengths_dict):
	cov_dict = {}
	for stat in stats:
		if len(stat) > 0 and \
			stat[0:5] != 'Total' and \
			stat.replace('\n', '').replace(' ', '') != '':
			name = stat.split('\t')[0]
			coverage = int(stat.split('\t')[1].replace('\n',''))
			if name not in lengths_dict.keys():
				continue
			coverage = coverage * READ_LEN/lengths_dict[name]
			# This /2 ensures that we are looking at averge 
			# read coverage in each replicate
			cov_dict[name] = coverage/2 
	return cov_dict

# Fill coverage dictionary with values from the RNAframework stats files
# RNAframework stats files is the number of reads that map onto each sequence
def fill_coverage_dict(stats_file, fasta_file):
	f = open(stats_file)
	stats = f.readlines()
	f.close()

	seq_lens = fill_length_dict(fasta_file)
	intron_cov_dict = get_norm_cov(stats, seq_lens)

	return intron_cov_dict

# Get reactivity values from xml files for an intron
def get_reac(cur_intron, reac_dir):
	cur_reac_file = os.path.join(reac_dir, cur_intron + '.xml')
	if os.path.exists(cur_reac_file):
		tree = ET.parse(cur_reac_file)
		reac_line = tree.getroot()[0][1].text.replace('\n', '')
		reac_line = reac_line.replace('\t', '')
		return reac_line.split(',')
	else:
		return []

def get_float_reac(cur_intron, reac_dir):
	reac = get_reac(cur_intron, reac_dir)
	return np.array([float(x) if x != 'NaN' else np.nan for x in reac])

def write_reactivity(name, reac_dir, cov, outdir, window_start, window_end):
	reacs = get_float_reac(name, reac_dir)

	if not os.path.exists(outdir):
		os.makedirs(outdir)

	if len(reacs) > 0:
		cur_tag = name.split("(")[0]
		f = open(outdir + '/' + name + '.txt', 'w')
		f.write('Coverage: %f\n' % cov[name])
		print("Coverage for %s: %f" % (name, cov[name]))
		reacs = reacs[window_start:(len(reacs) - window_end)]
		for reac in reacs:
			f.write('%f\n' % reac)
		f.close()

# This is for writing final reactivities that are an average of the two experiments
def write_all_reactivities(stats_file, fasta_file, reac_dir, outdir, window_start, window_end):
	names, seqs, name_seq_dict = get_names_seqs(fasta_file)
	cov_dict = fill_coverage_dict(stats_file, fasta_file)

	for cur_name in names:
		write_reactivity(cur_name, reac_dir, cov_dict, outdir, window_start, window_end)

write_all_reactivities(stats_file, fasta_file, reac_dir, outdir, window_start, window_end)

