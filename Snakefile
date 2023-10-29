ruleorder: Canu > minimap > Racon > tapestry > longstitch > teloclip > bwa_mem > sam2bam

indir = '/home/anouk/NOBINFBACKUP/Black_Sigatoka/data/'#'Cuban' #'data'

genomes = os.listdir(indir)
print(genomes)
genomes = 'P124'
"""
files_dict = {'C43':'Cuban/C43/C43.reads.fastq'}
genome_size_dict = {'C43':'46'}

"""#Fijiensis
files_dict = {g:f"{indir}/{g}/{os.listdir(f'{indir}/{g}')[0]}" for g in genomes}
#value = size in mbp (only number for flexibility)
genome_size_dict = {'P124': '42'}

print(files_dict)


rule all:
	input:
		expand('{genome}/{genome}_canu/{genome}.contigs.fasta', genome = genomes),
		expand('{genome}/{genome}_racon/{genome}_racon.fasta', genome = genomes),
		expand('{genome}/assembly1_{genome}', genome = genomes),
		expand('{genome}/longstitch/{genome}.scaffolds.fa', genome = genomes),
		expand('{genome}/teloclip_{genome}/{genome}vsreads.bam', genome = genomes),
		expand('{genome}/BWA/{genome}_illumina.sam', genome = genomes),
		expand('{genome}/BWA/{genome}_illumina.sorted.bam', genome = genomes),
		expand('{genome}/pilon/{genome}_pilon.fasta', genome = genomes)

rule Canu:
	input:
		lambda wildcards: files_dict[wildcards.genome]
	output:
		contigs_out = '{genome}/{genome}_canu/{genome}.contigs.fasta'
	params:
		outdir = '{genome}/{genome}_canu',
		genome_size = lambda wildcards : genome_size_dict[wildcards.genome]
	shell:
		'echo ~/anouk2/program/canu-2.2/build/bin/canu -p {wildcards.genome} \
		-d {params.outdir} genomeSize={params.genome_size}m -nanopore {input} \
		-OvlMerThreshold=300 -corMaxEvidenceErate=0.15'

rule minimap:
	input: 
		canudir = '{genome}/{genome}_canu/{genome}.contigs.fasta',
		rawread = lambda wildcards: files_dict[wildcards.genome]
	output:
		'{genome}/{genome}_racon/{genome}vsreads.sam'
	params:
		'{genome}/{genome}_canu/{genome}.contigs.fasta',
	shell:
		'minimap2 -ax map-ont {params} {input.rawread} > {output}'

rule Racon:
	input:
		sam = '{genome}/{genome}_racon/{genome}vsreads.sam',
		rawread = lambda wildcards: files_dict[wildcards.genome]
	output:
		'{genome}/{genome}_racon/{genome}_racon.fasta'
	params:
		'{genome}/{genome}_canu/{genome}.contigs.fasta',
	shell:
		'racon {input.rawread} {input.sam} {params} -t 30 > {output}'

rule tapestry:
	input:
		assembly = '{genome}/{genome}_canu/{genome}.contigs.fasta',
		reads = lambda wildcards: files_dict[wildcards.genome]
	output:
		'{genome}/assembly1_{genome}'
	shell:
		"weave -a {input.assembly} -r {input.reads} -t TTAGGG -o {output} -c 8"

rule longstitch:
	input:
		racon = '{genome}/{genome}_racon/{genome}_racon.fasta',
		reads = lambda wildcards: files_dict[wildcards.genome]
	output:
		link_racon = '{genome}/longstitch/{genome}.fa',
		link_reads = '{genome}/longstitch/{genome}.fq.gz',
		scaffold = '{genome}/longstitch/{genome}.scaffolds.fa'
	params:
		genome_size = lambda wildcards : genome_size_dict[wildcards.genome],
		directory = '{genome}/longstitch/',
	shell:
		"""
		mkdir -p {params.directory}
		ln -s /linuxhome/tmp/anouk/{input.racon} {output.link_racon}
		ln -s /linuxhome/tmp/anouk/{input.reads} {output.link_reads} 
		longstitch run -C {params.directory} out_prefix={wildcards.genome} \
		draft={wildcards.genome} reads={wildcards.genome} G={params.genome_size}e6
		"""

rule teloclip:
	input:
		scaffold = '/linuxhome/tmp/anouk/{genome}/longstitch/{genome}.scaffolds.fa',
		reads = lambda wildcards: files_dict[wildcards.genome]
	output:
		'{genome}/teloclip_{genome}/{genome}vsreads.bam'
	shell:
		"""
		samtools faidx {input.scaffold}
		minimap2 -ax map-ont {input.scaffold} {input.reads} \
		| samtools view -h -F 0x2308 \
		| teloclip --ref {input.scaffold}.fai | samtools sort > {output}
		"""

illumina_data = os.listdir('/media/anouk/KG000395/KG000395/fastq/')

rule bwa_index:
	input:
		assembly = '{genome}/teloclip_{genome}/{genome}_telo_compl.fasta',
	output:
		'{input.assembly}.fai'
	shell:
		'bwa index {input.assembly}'

rule bwa_mem:
	input: 
		assembly = '{genome}/teloclip_{genome}/{genome}_telo_compl.fasta',
		reads_1 = lambda wildcards: \
		[f'/media/anouk/KG000395/KG000395/fastq/{x}' for x in illumina_data if wildcards.genome in x and '_R1_' in x],
		reads_2 = lambda wildcards: \
		[f'/media/anouk/KG000395/KG000395/fastq/{x}' for x in illumina_data if wildcards.genome in x and '_R2_' in x],
	output:
		'{genome}/BWA/{genome}_illumina.sam'
	shell:
		"""
		bwa index {input.assembly} 
		bwa mem -t 10 -x ont2d {input.assembly} {input.reads_1} {input.reads_2} > {output}
		"""

rule sam2bam:
	input:
		'{genome}/BWA/{genome}_illumina.sam'
	output:
		bam = '{genome}/BWA/{genome}_illumina.sorted.bam',
		bam_index = '{genome}/BWA/{genome}_illumina.sorted.bam.bai'
	shell:
		"""
		samtools sort {input} -o {output.bam}
		samtools index {output.bam}
		"""

#rule bam_index:
#	input:
#		'{genome}/BWA/{genome}_illumina.sorted.bam'
#	output:
#		'{genome}/BWA/{genome}_illumina.sorted.bam.bai'
#	shell:
#		"samtools index {input} > {output}"

rule pilon:
	input:
		assembly = '{genome}/teloclip_{genome}/{genome}_telo_compl.fasta',
		bam = '{genome}/BWA/{genome}_illumina.sorted.bam',
	output:
		full = '{genome}/pilon/{genome}_pilon.fasta'
	params:
		directory = '{genome}/pilon',
		name = '{genome}_pilon',
	shell:
		"""
		pilon --genome {input.assembly} --bam {input.bam} --output {params.name} --outdir {params.directory} --changes --vcf --tracks --fix bases -Xmx200G 
		"""
