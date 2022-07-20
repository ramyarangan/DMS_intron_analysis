Produces Fig.2D-F and Fig. S4A-C.

run_m2seq_plots is a script to produce plots. Requires Biers package (https://ribokit.github.io/Biers/). The trim sequence ensures that reference hairpins are not included in the structure prediction.

Raw .rdat files from M2seq data analysis are provided.

Example usage: run_m2seq_plots( 'rdat_files/RTB008.reactivity.rdat', 'rdat_files/RTB018.reactivity.rdat', 'QCR9', 'figures/', 500, 0, 0, 27, 44, '')

