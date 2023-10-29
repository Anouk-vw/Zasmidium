"""
Script to get ITS from genome assembly and build the tree
"""

#1. Get ITS from assembly
ITSx -i ${genome_assembly} -o ${genome} -cpu 4 -t F

#2. Download ITS from NCBI > full fasta

#3. Align ITS_fasta mafft
mafft ${ITS_fasta} > ${alignment}

#4. Build tree from fasta
raxmlHPC -m GTRGAMMAI -p 5 -s ${alignment} -n ITS_tree_branchlength_full -N 500 -T 4 -k -f a -x 5
