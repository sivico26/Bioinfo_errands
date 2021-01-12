#!/usr/bin/env python
from time import time
from Bio import SeqIO
import argparse

def get_arguments():
	parser = argparse.ArgumentParser(description='Calculate number of reads, amount of information (bp), and mean read length for a folder with fasta/fastq',
			formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument('input', type=str,
			help='Path to file input')
	parser.add_argument('-f', '--format', type=str, default="fastq", choices=["fasta","fastq"],
			help='Input file format [%(default)s]')
#	parser.add_argument('-u', '--unique',action="store_true",
#                       help='Input is a single file [False]')
	args = parser.parse_args()
	return args

def get_read_stats(file,format):
	handle = SeqIO.parse(file,format)
	nreads = 0
	length = 0
	for read in handle:
		length += len(read)
		nreads += 1
	return length,nreads

def main():
	start_time = time()
	args = get_arguments() 
	length,nreads = get_read_stats(args.input,args.format)
	mean = length/nreads

	print('Amount of information: {} Mbp'.format(length/1000000))
	print('Number of reads: {}'.format(nreads))
	print('Mean read-length: {} bp'.format(mean))
	print('Execution time: {} seconds'.format(time() - start_time))
	
if __name__ == '__main__':
        main()



