#!/usr/bin/env python3
import argparse
import sys
from Bio import SeqIO

def parse_args():	
	parser = argparse.ArgumentParser(prog="filter_contigs",
		description=f"Filter short or bad contigs.")

	parser.add_argument('-i','--input', required=True,
		help='Input file in fasta format.')
	parser.add_argument('-l','--min-length', default=250,
		help='Minimum contig length. Default: 250.',type=int)
	parser.add_argument('-o','--output',
		help='Output file. If none print to stdout.')
	
	args = parser.parse_args()
	
	return args


def main():
	args = parse_args()
	i = 0
	if args.output:
		with open(args.output, 'w') as oh:
			for rec in SeqIO.parse(args.input, 'fasta'):
				if len(rec) >= args.min_length:
					SeqIO.write(rec, oh, 'fasta')
					i +=1
	else:
		for rec in SeqIO.parse(args.input, 'fasta'):
			if len(rec) >= args.min_length:
				SeqIO.write(rec, sys.stdout, 'fasta')
				i +=1
	print(f"{i} contigs passed filter", file=sys.stderr)

if __name__ == '__main__':
	main()