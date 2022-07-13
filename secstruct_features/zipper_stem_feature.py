import subprocess # For making RNAfold/RNAcofold calls
import os
import scipy.stats as stats

import matplotlib.pyplot as plt
from secstruct_util import *

# This code finds stable stems with dG < DG_LIM that fall within a region that is reasonable
# in the context of the spliceosome. The viable region must be pre-specified based 
# on structural analysis.

# Approach: Starts with reasonably long stems present in the MFE, and then checks the 
# stability of these stems and a slightly expanded surrounding context with RNAfold or RNAcofold. 
DG_LIM=0

# Region is (lowest start, highest end, total distance max, total distance min)
# Stem must start after lowest start + start_idx and must end before end_idx - highest end
FIRST_STEP_REGION = (10, 20, 12+10+20, 55+10+20) # For 5'SS to BP
SECOND_STEP_REGION = (0, 0, 0, 3000) # (8, 3, 8+3, 3000) # For BP to 3'SS

# Gets dictionary of start of BP: end of BP
def get_base_pairs(secstruct):
	ii = 0
	stack = []
	bps = {}
	while (ii < len(secstruct)):
		if (secstruct[ii] == "("):
			stack += [ii]
		if (secstruct[ii] == ")"):
			bps[stack[len(stack)-1]] = ii
			stack = stack[:(len(stack)-1)]
		ii += 1
	return bps

# Get stems of the form [(a, b, c, d, n)] for a---b paired to c---d 
# with n bps present in the stem, with max bulge size specified.
# Can handle strings with & in the middle to indicate two strands
def get_stems(secstruct, bulge_size=5):
	bps = get_base_pairs(secstruct)

	stems = []
	ii = 0
	# Search for the start of a new stem
	while (ii < len(secstruct)) and (secstruct[ii] != '&'):
		if ii not in bps:
			ii += 1
			continue
		start = ii
		end = bps[ii]
		jj = bps[ii]
		ii += 1
		bulge_cnt = 0
		n_bps = 1
		# Extend stem as far as possible, as long as no bulge is > limit
		while (ii < jj) and (bulge_cnt < bulge_size + 1):
			# End if you encounter a )
			if secstruct[ii] == ')' or secstruct[ii] == '&':
				break
			# Increment bulge size if necessary
			if secstruct[ii] == '.':
				bulge_cnt += 1
			# If a match is found, restart bulge count and increase bp count
			elif ((bps[ii] <= jj) and (bps[ii] >= jj - bulge_size)):
				jj = bps[ii]
				n_bps += 1
				# Set bulge_cnt = 0 to allow each junction to have bulge_size nts
				bulge_cnt = 0 
			else: # The endpoint is a different stem
				break

			# Next position
			ii += 1
		stems += [(start, ii - bulge_cnt, jj, end + 1, n_bps)]
	return stems

# Gets location of longest stem with bulge size limit of 3 within this secstruct
def get_max_stem(secstruct, bulge_size=3):
	stems = get_stems(secstruct, bulge_size=bulge_size)
	max_stem = (-1, -1, -1, -1, -1)
	for stem in stems:
		if (stem[4] > max_stem[4]):
			max_stem = stem
	return max_stem

# From the stem locations in 'stems', find all that are compatible with the 
# region limits specified in 'region'
# 'Region' is (lowest start, highest end, total distance max, total distance min)
# Stem must start after lowest start + start_idx and must end before end_idx - highest end
def get_stems_region(stems, region, start_idx, end_idx):
	matching_stems = []
	(low, high, max_dist, min_dist) = region
	#print(start_idx)
	#print(end_idx)
	for stem in stems:
		(a, b, c, d, n_bps) = stem
		#print("Stem constraints")
		#print(a - start_idx)
		#print(end_idx - d)
		if (a - start_idx) > low and (end_idx - d) > high and \
			((a - start_idx) + (end_idx - d)) > max_dist and \
			((a - start_idx) + (end_idx - d)) < min_dist:
			matching_stems += [stem]
	return matching_stems

def collect_dG_secstruct(seq, sys_command):
	f = open('tmp.dat', 'w')
	f.write(seq.replace('T', 'U'))
	f.close()
	try: 
		p = subprocess.Popen(sys_command + ' tmp.dat', shell=True, stdout=subprocess.PIPE)
		lines = p.stdout.readlines()
		os.remove('tmp.dat')

		# String parsing to process RNAfold/cofold output
		secstruct = lines[1].decode("utf-8").split()[0]
		# Handles cases when the output is like ((((..(((((...&)))))..)))) ( -8.20)\n or like
		# ((((..(((((...&)))))..)))) (-8.20)\n
		dG_str = ''.join(lines[1].decode("utf-8").split()[1:]) 
		dG = float(dG_str[1:-1])
	except: 
		dG = DG_LIM
		secstruct = ""
	return [dG, secstruct]

# System call to run RNAfold
def run_rnafold(seq):
	return collect_dG_secstruct(seq, 'RNAfold')

# System call to run RNAcofold
def run_cofold(seq):
	return collect_dG_secstruct(seq, 'RNAcofold')

# For an expanded region around the stem candidate from the MFE, check the stability
# using RNAfold or RNAcofold
def get_stem_dG(stem, seq, min_start, max_end, expand_low=5):
	
	# Get region surrounding the stem that does not go too close to the 5'SS and BP
	start1 = max(min_start, stem[0] - expand_low)
	end1 = stem[1] + expand_low
	end2 = min(max_end, stem[3] + expand_low)
	start2 = stem[2] - expand_low

	# Assemble input for RNAfold or RNAcofold
	run_seq = ""
	secstruct = ""
	dG = DG_LIM
	seq1 = seq[start1:min(end1, start2)]
	seq2 = seq[max(end1, start2):end2]
	run_seq = seq1 + "&" + seq2
	if len(run_seq) > 2:
		[dG, secstruct] = run_cofold(run_seq)

	stem_str = run_seq + "\n" + secstruct + "\n" + str(dG)
	return [dG, stem_str]

# Get the highest average BPP for a helix of length at least 4 within the stem
def get_best_bpp(secstruct, stem, bpp_matrix, len_cutoff=4, bpp_thresh=0.7):
	
	strand1 = secstruct[stem[0]:stem[1]]
	strand2 = secstruct[stem[2]:stem[3]]

	bps = get_base_pairs(secstruct)

	ii = stem[0]
	cur_end = stem[3] - 1
	cur_cnt = 0
	total_bpp = 0
	bpps = []
	cnts = []
	while ii < stem[1] and cur_end >= stem[2]:
		if (ii in bps) and (bps[ii] == cur_end):
			cur_cnt += 1
			# For now get the maximum base-pairing probability in the stem 
			# This is what is done to get the VARNA diagrams from Biers
			total_bpp = max(bpp_matrix[ii][cur_end], total_bpp)
			# print("%d %d %d %f\n" % (cur_cnt, ii, cur_end, bpp_matrix[ii][cur_end]))
			ii += 1
			cur_end -= 1
			continue

		if secstruct[ii] != '(':
			ii += 1
		elif cur_end != bps[ii]:
			cur_end -= 1

		if cur_cnt > 0:
			bpps += [total_bpp]
			cnts += [cur_cnt]

		cur_cnt = 0
		total_bpp = 0
	
	bpps += [total_bpp]
	cnts += [cur_cnt]

	#print(cnts)
	#print(bpps)
	# How many base-pairs are in stems passing the base-pair probability threshold?
	best_bpp = 0
	for ii, cnt in enumerate(cnts):
		#print(bpps[ii])
		if bpps[ii] > bpp_thresh:
			if cnt > len_cutoff: 
				best_bpp = max(bpps[ii], best_bpp)

	if best_bpp > 0:
		return True, best_bpp
	# total_bps_passing = 0
	# total_bpp = 0
	# for ii, cnt in enumerate(cnts):
	# 	#print(bpps[ii])
	# 	if bpps[ii] > bpp_thresh:
	# 		total_bps_passing += cnt
	# 		total_bpp = max(total_bpp, bpps[ii])
# 
	# if total_bps_passing > len_cutoff:
	# 	return True, total_bpp # /total_bps_passing

	# total_bps_passing = 0
	# total_bpp = 0
	# for ii, cnt in enumerate(cnts):
	#	print(bpps[ii]/cnt)
	#	if bpps[ii]/cnt > bpp_thresh:
	#		total_bps_passing += cnt
	#		total_bpp += bpps[ii]

	#if total_bps_passing > len_cutoff:
	#	return True, total_bpp/total_bps_passing
	return False, -1


def get_best_zipper_stem(bp, seq, mfe, bpp_matrix, min_num_bp=8, bpp_cutoff=0.7, \
	do_second=False, overhang_5p=0, overhang_3p=0):
	stems = get_stems(mfe)
	#print(stems)
	stems_region = []
	if do_second:
		stems_region = get_stems_region(stems, SECOND_STEP_REGION, bp, \
			len(seq) - overhang_3p) # Second step
	else:
		stems_region = get_stems_region(stems, FIRST_STEP_REGION, overhang_5p, bp) # First step
	#print(stems_region)

	has_zipper_stem = False
	best_dG = 200 # Some large number
	best_stem = ""
	for stem in stems_region:
		if stem[4] < min_num_bp:
			continue
		passes_bpp, best_bpp = get_best_bpp(mfe, stem, bpp_matrix)
		#print(best_bpp)
		if not passes_bpp:
			continue
		dG = 200
		cur_stem = ""
		if do_second:
			[dG, cur_stem] = get_stem_dG(stem, seq, bp + SECOND_STEP_REGION[0], \
				len(seq) - overhang_3p - SECOND_STEP_REGION[1]) # Second step
		else:
			[dG, cur_stem] = get_stem_dG(stem, seq, overhang_5p + \
				FIRST_STEP_REGION[0], bp - FIRST_STEP_REGION[1]) # First step
		if (dG < best_dG):
			has_zipper_stem = True
			best_dG = dG
			best_stem = stem

	return [has_zipper_stem, best_stem, best_dG]

def get_best_zipper_stem_tag(name, fasta_file, secstruct_dir, dat_file, do_second=False):
	_, _, name_seq_dict, _ = get_names_seqs_from_fasta(fasta_file)
	mfe = get_dotbracket(name, secstruct_dir)
	bpp_matrix = get_bpp_matrix(name, secstruct_dir)
	bp = get_bp_loc(name, dat_file)
	#print(mfe is not None)
	#print(bpp_matrix is not None)
	#print(bp)

	has_zipper_stem = False
	if (mfe is not None) and (bpp_matrix is not None) and (bp > 0):
		if do_second:
			has_zipper_stem, stem, best_dG = get_best_zipper_stem(bp, name_seq_dict[name], \
				mfe, bpp_matrix, min_num_bp=6, do_second=True)
		else:
			has_zipper_stem, stem, best_dG = get_best_zipper_stem(bp, name_seq_dict[name], \
				mfe, bpp_matrix, min_num_bp=6)

	chr_num = name.split(':')[0]
	chr_start = int(name.split(':')[1].split('-')[0])
	
	chr_num_str = ""
	strand1_start = ""
	strand1_end = ""
	strand2_start = ""
	strand2_end = ""
	dG = ""
	if has_zipper_stem:
		chr_num_str = chr_num
		strand1_start = str(chr_start + stem[0])
		strand1_end = str(chr_start + stem[1] - 1)
		strand2_start = str(chr_start + stem[2])
		strand2_end = str(chr_start + stem[3] - 1)
		dG = str(best_dG)

	#print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (name, chr_num_str, strand1_start, \
	#	strand1_end, strand2_start, strand2_end, dG))
	return (name, chr_num_str, strand1_start, strand1_end, strand2_start, strand2_end, dG)


def test_stems():
	#print(get_stems("...(((((..(((((....))......(((...)))..)))....)))))"))
	#name = 'chrIV:491559-491898'
	#_, _, name_seq_dict, _ = get_names_seqs_from_fasta("../intron_annot/standard_introns.fa")
	#mfe = get_dotbracket(name, "../reactivity/rnastructure_sherlock_1221/intron/")
	#bpp_matrix = get_bpp_matrix(name, "../reactivity/rnastructure_sherlock_1221/intron/")
	#print(get_best_zipper_stem(315, name_seq_dict[name], mfe, bpp_matrix))


	fasta_file = "../intron_annot/standard_introns.fa"
	secstruct_dir = "../reactivity/rnastructure_sherlock_1221/intron/"
	dat_file = "../intron_annot/standard_introns_allsize.dat"

	name = "chrXVI:654167-654570"
	print(get_best_zipper_stem_tag(name, fasta_file, secstruct_dir, dat_file, do_second=False))

