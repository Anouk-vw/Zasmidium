Repository containing the scripts used for assembly and analysis of the 
<i>Zasmidium syzygii</i> genome.

- Snakefile: 

	Script to create the genome assembly from Nanopore long reads. Prepares 
	all files for teleoclip to manually inspect missing telomeres.
- annotate.sh:

	Script to annotate the genome assembly
- quality_check.sh:

	Comments to check the quality of the genome assembly.
	(Genome size estimation, Genome completeness, Assembly statistics)
- ITS_analyis.sh:

	Comments to extract ITS from assembly and build phylogenetic tree.
