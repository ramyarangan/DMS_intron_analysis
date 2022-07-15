"""
Utility functions for getting sequences, coverage, secondary structures info,
reactivities, and splicing levels.
"""
import os
import xml.etree.ElementTree as ET
import numpy as np

snorna_introns = ['chrI:142253-142619', 'chrVII:364964-365432', 'chrXIII:499878-500151', 'chrXI:283095-283421',\
	'chrXII:856574-857057', 'chrXIII:163308-163716', 'chrXIV:721770-722302']

# Get the sequence tags and sequences from a fasta file in 
# a couple data structures (lists and a dictionary)
def get_names_seqs_from_fasta(fasta_file):
	f = open(fasta_file)
	lines = f.readlines()
	f.close()

	names = []
	seqs = []
	name_seq_dict = {}
	name_symbol_dict = {}
	cur_seq = ''
	for line in lines:
		if line[0] == '>':
			if cur_seq != '':
				seqs += [cur_seq]
				name_seq_dict[names[-1]] = cur_seq
			cur_seq = ''
			name = line.split(' ')[0].replace('>', '').replace('\n', '')
			gene_symbol = ''
			if len(line.split(' ')) > 1:
				gene_symbol = line.split(' ')[1].replace('\n', '')
			names += [name]
			name_symbol_dict[name] = gene_symbol
		else:
			cur_seq += line.replace('\n', '')
	seqs += [cur_seq]
	name_seq_dict[names[-1]] = cur_seq
	return names, seqs, name_seq_dict, name_symbol_dict

# Get coverage for an intron based on a tag (like chrI:11351-11535)
# Returns coverage (in RPKM)
def get_cov(name, reac_dir):
	reac_filename = reac_dir + name + ".txt"

	if not os.path.exists(reac_filename):
		return -1

	f = open(reac_filename)
	lines = f.readlines()
	f.close()

	return float(lines[0].split(' ')[1])

def get_dotbracket(name, secstruct_dir):
	dotbracket_filename = secstruct_dir + name + '_secstruct.txt'
	if not os.path.exists(dotbracket_filename):
		return None

	f = open(dotbracket_filename)
	lines = f.readlines()
	f.close()
	return lines[0].replace('\n', '')

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

def get_bp_loc(name, dat_file):
	# Standardize name tag format
	name = name.split("(")[0]
	
	f = open(dat_file)
	dat_lines = f.readlines()
	f.close()

	for ii in range(int(len(dat_lines)/2)):
		cur_bp = dat_lines[ii * 2].split(' ')[0]
		cur_items = dat_lines[ii * 2].split(' ')[1].split('\t')
		cur_tag = cur_items[0] + ":" + cur_items[1] + "-" + cur_items[2]
		if cur_tag == name:
			return int(cur_bp)

	return -1

def get_reactivities(name, reac_dir):
	reac_filename = reac_dir + name + ".txt"
	
	if not os.path.exists(reac_filename):
		return np.array([]), 0

	f = open(reac_filename)
	lines = f.readlines()
	f.close()

	cov = float(lines[0].split(' ')[1])
	reacs = []
	cnt = 0
	for line in lines[1:]:
		if line.replace('\n', '') == 'nan':
			reacs += [np.nan]
		else: 
			reacs += [float(line)]
		cnt += 1

	return reacs, cov

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

def get_ri_fractions(ri_fractions_file):
	f = open(ri_fractions_file)
	ri_fraction_lines = f.readlines()
	f.close()

	pladb_ri = {}
	no_pladb_ri = {}
	for line in ri_fraction_lines:
		tag = line.split(' ')[0][:-1]
		pladb_ri[tag] = float(line.split(' ')[1])
		no_pladb_ri[tag] = float(line.split(' ')[2])

	return pladb_ri, no_pladb_ri



