### **Fig. 2D-F, Fig. S4A-C, and Fig. S5**

run_m2seq_plots produces the M2-seq plots and secondary structures. Requires Biers package (https://ribokit.github.io/Biers/). The trim sequence positions ensure that reference hairpins are not included in the structure prediction, and should have the same values (27, 44) for all examples included here.

Raw .rdat files from M2seq data analysis are provided.

Example usage: run_m2seq_plots( 'rdat_files/RTB008.reactivity.rdat', 'rdat_files/RTB018.reactivity.rdat', 'QCR9', 'figures/', 500, 0, 0, 27, 44, '')

Example figures: 
[QCR9_Z_scores_plot.pdf](https://github.com/ramyarangan/DMS_intron_analysis/files/9147028/QCR9_Z_scores_plot.pdf)
![QCR9_VARNA](https://user-images.githubusercontent.com/2606810/179911899-d1469674-0e81-42d4-bca8-a739b8ca7421.png)
