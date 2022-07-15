# DMS_intron_analysis
Processing splicing inhibition and DMS-MaPseq data for introns and analyzing DMS-guided secondary structure predictions.

# Overview
This repository contains:
* Data from DMS-MaPseq after splicing inhibition:
    * Normalized reactivity data, raw mutational frequencies, and read coverage across different sequence sets in `analyze_run/combined_1221`
    * DMS-guided secondary structure predictions across different sequence sets in `secstruct_features/rnastructure_sherlock_1221/`
* Evaluating and processing DMS-MaPseq data in `analyze_run/`
    * Checking DMS modification patterns and DMS reactivity for rRNA
    * Evaluating the correlation between replicate experiments
    * Analyzing the impact of splicing inhibition by pladB treatment
    * Writing raw DMS-MaPseq output to convenient files 
* Evaluating DMS structure predictions for control structures in `control_structs/`
    * Control structure native structures in `control_structs/native_secstructs/`
    * Predicted structures in `control_structs/control_secstruct_bpps/`
    * Evaluating PPV and sensitivity for structure predictions with different approaches
* Sequence files for different intron / control sequence sets in `intron_annot/`
    * .fa for fasta format
    * .dat for sequence + branchpoint and strand direction)
* Secondary structure analysis for DMS-guided secondary structure predictions in `secstruct_features/`
    * Compute features: Gini index, zipper stems and end stems, longest stem, maximum path length, average DMS accessibility
    * Compare features between introns and coding regions: statistical comparisons and plots
    * Cluster introns by secondary structure features: hierarchical clustering, tSNE, PCA
    * Write table compiling intron information, secondary structures, and structural features

# Requirements
* xml, scipy, numpy, matplotlib, seaborn, pandas
* arnie, Vienna
