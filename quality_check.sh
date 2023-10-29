genome='P124'
genome_assembly='path/to/genomeassembly/'
Illumina_reads='path/to/reads'
short_read_bam='path/to/shortreadbam/'
blast_db='path/to/blastdb/'

#Predict genome size Jellyfish doi:10.1093/bioinformatics/btr011
jellyfish count -C -m 21 -s 1000000 -t 4 ${Illumina_reads} -o ${genome}.jf
jellyfish histo -t 4 ${genome}.jf > ${genome}.histo

#Plot short read coverage
java -jar /path/to/jvarkit/dist/jvarkit.jar wgscoverageplotter -R ${genome_assembly} ${short_read_bam} --min-contig-length 10kb > WGS_coverage/${genome}_wgs_plot.svg

#BLOBtools
#run blastn (long run time)
blastn -db ${blast_db} -query ${genome_assembly} -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads 16 -out ${genome}.blast.out

blobtools create --fasta ${genome_assembly} -b ${short_read_bam} -t blast 
