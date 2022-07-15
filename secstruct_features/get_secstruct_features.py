"""
Get secondary structure features for a given sequence tag including: 
* longest stem
* maximum path length
* average stem helix confidence estimates
* number of NWJs

Constructs two types of graph representations for a secondary structure:
- Fine-grained graph with base-pairs and single-stranded nucleotides as nodes
- Coarse-grained graph with stems and junctions as nodes
"""
import sys
import networkx as nx 
import os
from matplotlib import pyplot as plt 
from scipy import stats 
import pandas as pd 
import numpy as np
import seaborn as sns 

from secstruct_util import *

COVERAGE_CUTOFF = 1971

class Node():
	def __init__(self, base1, base2=-1, base_pair=False):
		self.base1 = base1
		self.base2 = base2
		self.base_pair = base_pair

class BigNode():
	def get_type(self):
		return 'base'

class Stem(BigNode):
	def __init__(self, strand1_nts, strand2_nts, bpp_matrix):
		if len(strand1_nts) != len(strand2_nts):
			raise RuntimeError("Stem needs to have equal # nts in strands")

		self.strand1_nts = strand1_nts
		self.strand2_nts = strand2_nts
		self.bpp = self.get_bpp(bpp_matrix)

	def get_bpp(self, bpp_matrix):
		total_bpp = 0
		for ii, nt1 in enumerate(self.strand1_nts):
			nt2 = self.strand2_nts[ii]
			total_bpp += bpp_matrix[nt1][nt2]
		return total_bpp/len(self.strand1_nts)

	def get_type(self):
		return 'stem'

	def len(self):
		return len(self.strand1_nts)

	def get_min(self):
		return min(min(self.strand1_nts), min(self.strand2_nts))

	def get_max(self):
		return max(max(self.strand1_nts), max(self.strand2_nts))

	def __str__(self):
		print_str = "Stem of length %d containing base-pairs:\n" % len(self.strand2_nts)
		for ii, n1 in enumerate(self.strand1_nts):
			print_str += "%d %d\n" % (n1, self.strand2_nts[ii])
		print_str += "Base pair probability: %f\n" % self.bpp
		return print_str

class Junction(BigNode):
	def __init__(self, nts):
		self.nts = nts

	def get_min(self):
		return min(self.nts)

	def get_max(self):
		return max(self.nts)

	def get_type(self):
		return 'junction'

class External(BigNode):
	def __init__(self, nts):
		self.nts = nts

	def get_min(self):
		return min(self.nts)

	def get_max(self):
		return max(self.nts)

	def get_type(self):
		return 'external'

# Node types: 
#    Base-pair (base-pairing residues)
#    Single nucleotide
# Used for: MLD calculation (longest shortest paths)
# Return type: networkx graph, node id to node object dict
def get_graph(dotbracket, ss_weight=0, bp_weight=1):
	base1_list = []
	node_dict = {}

	num_nodes = 0
	for ii, curchar in enumerate(dotbracket):
		if curchar == '.':
			num_nodes += 1
			node_dict[ii] = Node(ii)
		if curchar == '(':
			base1_list += [ii]
		if curchar == ')':
			base1 = base1_list[-1]
			base1_list = base1_list[:-1]
			new_node = Node(base1, base2=ii, base_pair=True)
			node_dict[ii] = new_node
			node_dict[base1] = new_node
			num_nodes += 1

	G = nx.Graph()

	for ii in range(len(dotbracket)):
		if node_dict[ii].base_pair:
			# Base pair node
			if node_dict[ii].base1 == ii:
				G.add_node(ii)
		else:
			# Loop node
			G.add_node(ii)
		
		if ii > 0:
			prev_node = node_dict[ii - 1]
			weight = ss_weight
			if node_dict[ii].base_pair and prev_node.base_pair:
				weight = bp_weight
			if not G.has_edge(prev_node.base1, node_dict[ii].base1):
				G.add_edge(prev_node.base1, node_dict[ii].base1, weight=weight)

	return G, node_dict

# Update graph with nodes for (cur_start, cur_end - 1) (inclusive)
# nt_to_node_dict: for each nucleotide, the stem, junction, or external node it belongs to
# dotbracket: dot bracket notation for the full structure
# G: current graph
# is_internal: False for a portion of the structure that is not in any stems 
#              or internal to any stems
def get_stem_graph_rec(cur_start, cur_end, nt_to_node_dict, \
	dotbracket, G, is_internal, verbose=False):
	# Add node for current portion of stem graph
	G.add_node(cur_start)

	if dotbracket[cur_start] == ')':
		raise RuntimeError("Unexpected recursion architecture: start index is a closing base-pair")

	# Case 1: stems. Here cur_start and cur_end are the ends of a stem; recursively update internal loops
	if dotbracket[cur_start] == '(':
		if verbose: 
			print("Doing stem at position: %d" % cur_start)
		if cur_start not in nt_to_node_dict.keys() or \
			nt_to_node_dict[cur_start].get_type() != 'stem':
			raise RuntimeError("Expected base-paired residue to be in a stem.")

		cur_stem = nt_to_node_dict[cur_start]
		
		stem_start = min(cur_stem.strand1_nts)
		junc_start = max(cur_stem.strand1_nts)
		junc_end = min(cur_stem.strand2_nts)
		stem_end = max(cur_stem.strand2_nts)
		
		if stem_start != cur_start:
			raise RuntimeError("Unexpected recursion architecture")

		# Update graph for nts between the 5' and 3' strands of stem
		# E.g. when Called on (((....)))... from (((((....)))...)), get node for ....
		node_id = get_stem_graph_rec(junc_start + 1, junc_end, nt_to_node_dict, \
			dotbracket, G, True, verbose=verbose)
		# Connect node for .... to ((()))
		G.add_edge(stem_start, node_id)
		
		node_id_2 = cur_start
		if cur_end > stem_end + 1:
			# This happens for internal loops with only 3' nucleotides or
			# external ssRNA outside of stems 
			# E.g. when Called on (((....)))... from (((((....)))...)), get node for ...
			node_id_2 = get_stem_graph_rec(stem_end + 1, cur_end, nt_to_node_dict, \
				dotbracket, G, is_internal, verbose=verbose)
			# Connect node for ... to ((()))
			G.add_edge(stem_start, node_id_2)

		if verbose:
			print("Finished stem at position: %d" % cur_start)

		if is_internal:
			# Called on (((....)))... from (((((....)))...))
			# Should return node for ... 
			# Called on (((....))) from ((...(((....)))...)) 
			# Should return node for ((()))
			return node_id_2
		# Called on (((....)))..... from ...(((....))).....
		# Should return node for ((()))
		return cur_start

	# all ssRNA in the loop / junctions / external nucleotides
	junc_nts = []

	# Case 2: external ssRNA
	if not is_internal:
		if verbose:
			print("Doing external ssRNA at position: %d" % cur_start)
		ssrna_end = cur_end
		cur_dotbracket = dotbracket[cur_start:cur_end]
		if '(' in cur_dotbracket: 
			ssrna_end = cur_dotbracket.index('(') + cur_start
			if ')' not in cur_dotbracket: 
				raise RuntimeError("Unexpected stem architecture")
			if dotbracket.index(')') < dotbracket.index('('):
				raise RuntimeError("Unexpected stem architectures")

		junc_nts = list(range(cur_start, ssrna_end))

		if cur_end > ssrna_end:
			# E.g. when called on ....(((...))) from ((....))....(((...)))
			node_id = get_stem_graph_rec(ssrna_end, cur_end, nt_to_node_dict, \
				dotbracket, G, False, verbose=verbose)
			# Connect .... to ((()))
			G.add_edge(cur_start, node_id)

		if verbose:
			print("Finished external ssRNA at position: %d" % cur_start)

		# Add new junction to datastructure
		new_external = External(junc_nts)
		for ii in junc_nts:
			nt_to_node_dict[ii] = new_external	

	# Case 3: internal ssRNA - loops and junctions
	else:
		if verbose:
			print("Doing internal loop / junction at position: %d" % cur_start)

		junc_nts = [] 
		neighbor_stems = [] # all stems internal to the junction
		cur_idx = cur_start
		while (cur_idx < cur_end):
			if dotbracket[cur_idx] == '.':
				junc_nts += [cur_idx]
				cur_idx += 1
			else:
				if cur_idx not in nt_to_node_dict.keys() or \
					nt_to_node_dict[cur_idx].get_type() != 'stem':
					raise RuntimeError("Expected base-paired residue to be in a stem.")
				neighor_stem = nt_to_node_dict[cur_idx]
				neighbor_stems += [neighor_stem]
				
				# Skip over all nucleotides in the new stem
				cur_idx = max(neighor_stem.strand2_nts) + 1

		# Recursion over all stems in the junction
		for neighbor_stem in neighbor_stems:
			start_idx = min(neighbor_stem.strand1_nts)
			end_idx = max(neighbor_stem.strand2_nts)
			# E.g. when called on ....((...))...(....).. from ((....((...))...(....)..))
			node_id = get_stem_graph_rec(start_idx, end_idx + 1, nt_to_node_dict, \
				dotbracket, G, True, verbose=verbose)
			# Connect ......... to (()) and ......... to ()
			G.add_edge(cur_start, node_id)

		if verbose:
			print("Finished internal loop / junction at position: %d" % cur_start)

		# Add new junction to datastructure
		new_junction = Junction(junc_nts)
		for ii in junc_nts:
			nt_to_node_dict[ii] = new_junction

	# When called on ....((...))...(....).. from ((....((...))...(....)..)) 
	# return start of .........
	# When called on ....(((...))) from ((....))....(((...)))
	# return start of ....
	return cur_start

# Node types: 
#      Stem (residues, length, bootstrapping probability) 
#      Loop (all loop nts condensed into 1)
# Used for: NWJ accounting, longest stem
# Return type: networkx graph, dictionary connecting nucleotides to the stems/ junctions
#              they belong in 
def get_stem_graph(dotbracket, bpp_matrix, stem_verbose=False, graph_verbose=False):
	stems = []

	# First build list of stems
	base1_list = []
	cur_stem = []
	nt_to_stem_dict = {}
	for ii, curchar in enumerate(dotbracket):
		in_stem = False
		if len(cur_stem) > 0 and curchar == ')':
			# We could be continuing a stem
			if base1_list[-1] == cur_stem[-1][0] - 1:
				in_stem = True

		if not in_stem and len(cur_stem) > 0:
			strand1_nts = [x[0] for x in cur_stem]
			strand2_nts = [x[1] for x in cur_stem]
			new_stem = Stem(strand1_nts, strand2_nts, bpp_matrix)
			if stem_verbose:
				print(new_stem)
			stems += [new_stem]
			for nt1 in strand1_nts:
				nt_to_stem_dict[nt1] = new_stem
			for nt2 in strand2_nts:
				nt_to_stem_dict[nt2] = new_stem
			cur_stem = []

		if curchar == '(':
			base1_list += [ii]

		if curchar == ')':
			cur_stem += [(base1_list[-1], ii)]
			base1_list = base1_list[:-1]

	if len(cur_stem) > 0:
		strand1_nts = [x[0] for x in cur_stem]
		strand2_nts = [x[1] for x in cur_stem]
		new_stem = Stem(strand1_nts, strand2_nts, bpp_matrix)
		stems += [new_stem]		
		for nt1 in strand1_nts:
			nt_to_stem_dict[nt1] = new_stem
		for nt2 in strand2_nts:
			nt_to_stem_dict[nt2] = new_stem

	# Initiate stem graph building recursion
	G = nx.Graph()
	nt_to_node_dict = nt_to_stem_dict
	get_stem_graph_rec(0, len(dotbracket), nt_to_node_dict, dotbracket, G, \
		False, verbose=graph_verbose)
	return G, nt_to_node_dict

def get_longest_stem(stemG, nt_to_node_dict, seq_len, bpp_cutoff=0.7, \
	loop_cutoff=10, bpp_len_cutoff=4, overhang_5p=0, \
	overhang_3p=0):

	# All junctions must be 2WJ
	# < 10 loop nucleotides
	tested_stems = set()
	longest_stem = []
	longest_len = 0
	for nt in nx.nodes(stemG):
		# Only explore each node once
		if nt in tested_stems: 
			continue
		tested_stems.add(nt)

		if nt not in nt_to_node_dict.keys():
			raise RuntimeError("Unexpected missing BigNode for nt: %s" % nt)

		# Extend current stem as far as possible
		cur_longest_stem = []
		cur_longest_len = 0

		if nt_to_node_dict[nt].get_type() != 'stem':
			continue
		
		cur_stem = nt_to_node_dict[nt]

		if cur_stem.get_min() < overhang_5p or \
			cur_stem.get_max() > seq_len - overhang_3p:
			continue

		if cur_stem.bpp < bpp_cutoff:
			continue

		passes_bpp_len = False
		if cur_stem.len() > bpp_len_cutoff:
			passes_bpp_len = True
		
		cur_longest_stem += [cur_stem]
		cur_longest_len += cur_stem.len()

		# Extend stem as far as possible
		neighbors = list(nx.neighbors(stemG, nt))
		while len(neighbors) > 0:
			cur_neighbor = neighbors[0]
			neighbors = neighbors[1:]

			if cur_neighbor in tested_stems:
				continue
			
			tested_stems.add(cur_neighbor)

			neighbor_node = nt_to_node_dict[cur_neighbor]


			if neighbor_node.get_min() < overhang_5p or \
				neighbor_node.get_max() > seq_len - overhang_3p:
				continue

			if neighbor_node.get_type() == 'stem':
				if neighbor_node.bpp < bpp_cutoff:
					continue
				if neighbor_node.len() > bpp_len_cutoff:
					passes_bpp_len = True
				cur_longest_stem += [neighbor_node]
				cur_longest_len += neighbor_node.len()

			if neighbor_node.get_type() == 'junction':
				if len(neighbor_node.nts) > loop_cutoff:
					continue
				next_neighbors = list(nx.neighbors(stemG, cur_neighbor))
				if len(next_neighbors) > 2:
					continue

			if neighbor_node.get_type() == 'external':
				continue
			# Add all neighbors to the neighbor list to search for
			# extensions to the stem. 
			new_neighbors = nx.neighbors(stemG, cur_neighbor)
			for candidate in new_neighbors:
				if candidate not in tested_stems:
					neighbors += [candidate]

		if passes_bpp_len and cur_longest_len > longest_len:
			longest_len = cur_longest_len
			longest_stem = cur_longest_stem

	return longest_stem, longest_len

def get_nwj(stemG, nt_to_node_dict, bpp_cutoff=0.9, \
	bpp_cutoff_2=0.5, loop_cutoff=20, min_n=3):
	nwjs = []

	for nt in nx.nodes(stemG):
		if nt not in nt_to_node_dict.keys():
			raise RuntimeError("Unexpected missing BigNode for nt: %s" % nt)

		if nt_to_node_dict[nt].get_type() != 'junction':
			continue

		# Check that there are at least min_n neighbors
		if len(list(nx.neighbors(stemG, nt))) != min_n:
			continue

		# Check that loop is at most loop_cutoff length
		if len(nt_to_node_dict[nt].nts) > loop_cutoff:
			continue

		# Check that 0 neighboring stems violate bpp_cutoff_2, and
		# at least 2 pass bpp cutoff 1
		num_bpp_cutoff_pass = 0
		num_bpp_cutoff_2_fail = 0
		for neighbor in nx.neighbors(stemG, nt):
			if neighbor not in nt_to_node_dict.keys():
				raise RuntimeError("Unexpected missing BigNode for nt: %s" % neighbor)
			neighbor_node = nt_to_node_dict[neighbor]
			if neighbor_node.get_type() != 'stem':
				raise RuntimeError("Unexpected loop/external node neighboring loop")
			if neighbor_node.bpp > bpp_cutoff:
				num_bpp_cutoff_pass += 1
			if neighbor_node.bpp < bpp_cutoff_2:
				num_bpp_cutoff_2_fail += 1
		#print("Num passing bpp_cutoff: %d; num failing cutoff 2: %d" % \
		#	(num_bpp_cutoff_pass, num_bpp_cutoff_2_fail))

		if num_bpp_cutoff_2_fail > 0: 
			continue
		if num_bpp_cutoff_pass < 2: 
			continue

		nwjs += [nt_to_node_dict[nt]]

	return nwjs

def get_mpl_from_graph(G, node_dict, norm_len, overhang_5p=0, overhang_3p=0):
	path_lens = nx.floyd_warshall(G, weight='weight')
	max_value = 0
	for cur_key in path_lens.keys():
		if cur_key < overhang_5p or cur_key > norm_len - overhang_3p:
			continue
		cur_max = 0
		for cur_key2 in path_lens[cur_key].keys():
			if cur_key2 < overhang_5p or cur_key2 > norm_len - overhang_3p:
				continue
			cur_max = max(cur_max, path_lens[cur_key][cur_key2]) 
		max_value = max(max_value, cur_max)

	return max_value/(norm_len - overhang_5p - overhang_3p)

def get_bpps(name, secstruct_dir, len_cutoff=6, overhang_5p=0, \
	overhang_3p=0):
	dotbracket = get_dotbracket(name, secstruct_dir)
	if dotbracket is None:
		return None
	bpp_matrix = get_bpp_matrix(name, secstruct_dir)
	if bpp_matrix is None:
		return None
	stemG, nt_to_node_dict = get_stem_graph(dotbracket, bpp_matrix, \
		stem_verbose=False, graph_verbose=False)

	all_bpp = []
	seen_nts = set()
	for cur_nt in nt_to_node_dict.keys():
		if cur_nt in seen_nts:
			continue
		cur_node = nt_to_node_dict[cur_nt]
		if cur_node.get_type() != 'stem':
			continue
		if cur_node.get_min() < overhang_5p or \
			cur_node.get_max() > len(dotbracket) - overhang_3p:
			continue
		for nt in cur_node.strand1_nts:
			seen_nts.add(nt)
		for nt in cur_node.strand2_nts:
			seen_nts.add(nt)
		if len(cur_node.strand1_nts) >= len_cutoff:
			all_bpp += [cur_node.bpp]
	return all_bpp


def get_bpps_tag(tag, fasta_file, secstruct_dir, cov_reac_dir, len_cutoff=6):
	names, seqs, _, _ = get_names_seqs_from_fasta(fasta_file)

	for ii, name in enumerate(names):
		if name == tag:
			coverage = get_cov(name, cov_reac_dir)

			if coverage > COVERAGE_CUTOFF:
				bpps = get_bpps(name, secstruct_dir, len_cutoff=len_cutoff)
				return bpps

			else: 
				return None

	return all_bpps

def get_longest_stem_len(tag, fasta_file, secstruct_dir, cov_reac_dir):
	names, seqs, _, _ = get_names_seqs_from_fasta(fasta_file)
	coverage = get_cov(tag, cov_reac_dir)

	if coverage < COVERAGE_CUTOFF:
		return -1

	longest_stem_len = -1

	for name in names:
		if name == tag:
			dotbracket = get_dotbracket(name, secstruct_dir)
			if dotbracket is None:
				return -1			
			bpp_matrix = get_bpp_matrix(name, secstruct_dir)
			if bpp_matrix is None:
				return -1

			stemG, nt_to_node_dict = get_stem_graph(dotbracket, bpp_matrix, \
				stem_verbose=False, graph_verbose=False)			
			_, longest_stem_len = get_longest_stem(stemG, nt_to_node_dict, len(dotbracket))

	return longest_stem_len


def get_num_nwj(tag, fasta_file, secstruct_dir, cov_reac_dir):
	names, seqs, _, _ = get_names_seqs_from_fasta(fasta_file)
	coverage = get_cov(tag, cov_reac_dir)

	if coverage < COVERAGE_CUTOFF:
		return -1

	num_nwj = -1

	for name in names:
		if name == tag:
			dotbracket = get_dotbracket(name, secstruct_dir)
			if dotbracket is None:
				return -1			
			bpp_matrix = get_bpp_matrix(name, secstruct_dir)
			if bpp_matrix is None:
				return -1
			stemG, nt_to_node_dict = get_stem_graph(dotbracket, bpp_matrix, \
				stem_verbose=False, graph_verbose=False)			
			num_nwj = len(get_nwj(stemG, nt_to_node_dict))

	return num_nwj

def get_mpl(name, seq, secstruct_dir, cov_reac_dir, overhang_5p=0, overhang_3p=0):
	coverage = get_cov(name, cov_reac_dir)
	
	dotbracket = get_dotbracket(name, secstruct_dir)
	if dotbracket is None:
		return -1
	G, node_dict = get_graph(dotbracket, ss_weight=1, bp_weight=1)
	return get_mpl_from_graph(G, node_dict, len(dotbracket), overhang_5p=overhang_5p, \
		overhang_3p=overhang_3p)

def get_mpl_tag(tag, fasta_file, secstruct_dir, cov_reac_dir):
	names, seqs, _, _ = get_names_seqs_from_fasta(fasta_file)
	coverage = get_cov(tag, cov_reac_dir)

	if coverage < COVERAGE_CUTOFF:
		return -1

	mpl = -1
	for name in names:
		if name == tag:
			dotbracket = get_dotbracket(name, secstruct_dir)
			if dotbracket is None:
				return -1
			G, node_dict = get_graph(dotbracket, ss_weight=1, bp_weight=1)
			mpl = get_mpl_from_graph(G, node_dict, len(dotbracket))
			break

	return mpl

def get_mpl_from_file(tag, mpl_stats_file, cov_reac_dir):	
	f = open(mpl_stats_file)
	mpl_lines = f.readlines()
	f.close()

	mpl = -1
	
	coverage = get_cov(tag, cov_reac_dir)

	if coverage < COVERAGE_CUTOFF:
		return mpl

	for ii in range(int(len(mpl_lines)/2)):
		name = mpl_lines[2 * ii].replace('\n', '')
		if name == tag:
			mpl = float(mpl_lines[2 * ii + 1])
			break

	return mpl

def get_length(tag, fasta_file):
	_, _, name_to_seq_dict, _ = get_names_seqs_from_fasta(fasta_file)

	return len(name_to_seq_dict[tag])
