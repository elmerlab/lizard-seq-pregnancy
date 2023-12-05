#!/bin/bash
# code snippet for extracting a tsv of gene SYMBOLs and protein ID
# by Hongxin

# get protein ID for each gene ID
grep -oE "gene=[^\;]+.*\;protein_id=[^\;]+" ../02_reference_data/annotation.gff|awk -F '\;' '{print $1, $NF}'|sort|uniq|sed 's/gene=//g'|sed 's/ protein_id=/\t/g' > ../02_reference_data/gene2prot.txt

# end