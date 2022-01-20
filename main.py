import sys
import os
from subprocess import call
import argparse
from zipfile import ZipFile
import pandas as pd
from Bio import SeqIO
import pysam

def parse_zip(zip_name):
	"""Get fastq and xls filenames and extract them"""
	with ZipFile(args.genexus_archive) as zf:
		zipped_files = zf.namelist()
		if 'COVID19AnnotateSnpEff.zip' not in zipped_files:
			sys.exit("COVID19AnnotateSnpEff.zip is absent")

		fastq_name = ""
		for filename in zipped_files:
			if filename.endswith('.fastq'):
				fastq_name = filename

		if not fastq_name:
			sys.exit("fastq file is absent")

		with open(fastq_name, 'wb') as f_fastq:
			f_fastq.write(zf.read(fastq_name))

		with open('COVID19AnnotateSnpEff.zip', 'wb') as f_zip:
			f_zip.write(zf.read('COVID19AnnotateSnpEff.zip'))

	with ZipFile('COVID19AnnotateSnpEff.zip') as zf_SnpEff:
		zipped_files = zf_SnpEff.namelist()
		xls_filename = ""
		for filename in zipped_files:
			if filename.endswith('.xls'):
				xls_filename = filename

		if not xls_filename:
			sys.exit("xls file is absent")

		with open(xls_filename, 'wb') as f_xls:
			f_xls.write(zf_SnpEff.read(xls_filename))

	return fastq_name, xls_filename

def get_coverage(fastq_name, reference, reference_len, chrom_name, cores_num):
	"""Get number of reads that cover each position"""
	if call(f"bowtie2-build {reference} covid_ref", shell=True)!=0:
		sys.exit(f"bowtie2-build for {reference} failed")

	if call(f"bowtie2 -p {cores_num} -x covid_ref -U {fastq_name} -S aligned.sam", shell=True)!=0:
		sys.exit(f"bowtie2 for {fastq_name} failed")

	if call(f"samtools sort -@ {cores_num} aligned.sam > aligned.bam", shell=True)!=0:
		sys.exit(f"sam to bam failed")

	if call(f"samtools index aligned.bam", shell=True)!=0:
		sys.exit(f"samtools index failed")

	bam = pysam.AlignmentFile('aligned.bam', 'rb')
	coverage = bam.count_coverage(chrom_name, 0, reference_len, quality_threshold = args.quality_threshold)
	zip_coverage = list(zip(coverage[0], coverage[1], coverage[2], coverage[3]))
	return [sum(zip_coverage[pos]) for pos in range(reference_len)]

if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='main.py')
	parser.add_argument('--genexus_archive', help='Genexus output zip file', required=True)
	parser.add_argument('--reference', help='Reference genome fasta file', required=True)
	parser.add_argument('--cores_num', help='Cores num for bowtie2, default = 16', default=16, type=int)
	parser.add_argument('--min_coverage', help='Minimum coverage per base, default = 1', default=1, type=int)
	parser.add_argument('--quality_threshold', help = 'Minimum quality score (in phred) a base has to reach to be counted, default=15', default=15, type=int)
	parser.add_argument('--output', help = 'Output fasta file name', required=True)
	args = parser.parse_args()
	if not args.genexus_archive:
		sys.exit("Provide genexus archive filename")

	# parse archive
	fastq_name, xls_filename = parse_zip(args.genexus_archive)

	# exclude frameshift_variant rows from snp_table
	snp_table = pd.DataFrame()
	xls = pd.read_csv(xls_filename, sep = '\t')
	snp_table = xls[['POS', 'REF', 'ALT', 'VARTYPE', 'ANN[*].EFFECT']]
	snp_table = snp_table.loc[snp_table['ANN[*].EFFECT'] != 'frameshift_variant']

	# read reference fasta
	reference_fasta = list(SeqIO.parse(args.reference, "fasta"))
	reference_seq = reference_fasta[0].seq
	consensus = list(reference_seq)
	reference_len = len(reference_seq)
	chrom_name = reference_fasta[0].id

	# get coverage for each position
	counts = get_coverage(fastq_name, args.reference, reference_len, chrom_name, args.cores_num)

	# replace positions with zero reads mapped with N
	nucl_dict = {0:'A', 1:'C', 2:'G', 3:'T'}
	reference_len = len(counts)
	for pos in range(reference_len):
		reads_num = counts[pos]
		if reads_num < args.min_coverage:
			consensus[pos] = 'N'

	# apply ALT if there is no Ns
	for index, row in snp_table.iterrows():
		pos = row['POS']
		replace_length = len(row['ALT'])
		if 'N' not in consensus[pos:pos+replace_length]:
			consensus = consensus[:pos]+list(row['ALT'])+consensus[pos+replace_length:]

	# get fasta header from id in zip filename
	path, file = os.path.split(args.genexus_archive)
	fasta_header = file.split('_')[0]
	
	with open(args.output, 'w') as f_output:
		f_output.write(f">{fasta_header}\n{''.join(consensus)}\n")