# TKO_CRISPR_screen
Analysis of the TKO CRISPR screen

1. [Extract & match guides](count_guides.TKO.pl) - this script extracts 20nt long guide sequences from the R1 of sequencing reads, using matches to the U6 promoter and gRNA backbone sequences with the pattern:  m/CACCG(.{20})GTTTTAGAGC/; then finds exact matches to the sequences from TKOv3 library

2. [Count tables](guide_count_table.R) - this script finds the guide summary files generated with the Perl script \*.guides.txt in the current working directory. Then merges them together into gene count tables, and guide tables, and normalizes them to the total number of detected guides.

3. [Counts distribution & log2FC changes](TKOvs_analysis.R) - explore distribution of counts per gene and per guide; compare REF1 with TKO library, and HIGH vs LOW insulin samples, etc.
