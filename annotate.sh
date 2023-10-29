"""
Annotate genome using Z. cellare protein evidence
"""

#mask repeats
RepeatMasker -lib ${repeat_consensus} ${genome_assembly} -xsmall -dir .


#Funannotate 
funannotate predict -i ${masked_assembly} -o fun_P124 -s Zasmidium_cellare\
  --busco_db pezizomycotina --cpus 1 --protein_evidence ${Cellare_proteins}




