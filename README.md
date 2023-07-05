Ribo-seq metadata standardization for RiboCrypt and Riboseq.org
==============================================================================

The script does 3 things:

0. Given a Entrez fetch table ~ 700 columns from SRA (not included in the scripts)
1. Standardize column names (CELL_LINE, CELL LINE, celllines are all the same)
2. Standardize column values: (Ribo-seq, Riboseq, RIBOSEQ are all the same)
3. Semi manual annotation (HeLA is female cell line, HEK is male etc)

Finally upload this file to google drive with statistics of how much
could be standardized. 

#### About

The procedure is packaged into 3 scripts ran from: metadata_main_script.R

