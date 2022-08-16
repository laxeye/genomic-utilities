#!/usr/bin/env python
import argparse
import re
import sys
import os.path

from Bio import SeqIO


def parse_args():
	parser = argparse.ArgumentParser(description="Contig renamer. "
		+ "Example run: renamer.py -i contigs/B-2335.fna -t \"Bacillus velezensis VKPM=B-2335\"")
	parser.add_argument('-i', '--input', required=True,
		help='Input file')
	parser.add_argument('-c', '--locus-tag',
		help='Locus tag, by default input filename.')
	parser.add_argument('-o', '--output',
		help='Output filename, by default same as input.')
	parser.add_argument('-m', '--min-len', type=int, default=250,
		help='Min length of contig to keep it, default 250 bp.')
	parser.add_argument('-t', '--taxonomy', default='Unknown sp.',
		help='Taxonomy string, eg. \"Bacillus sp. VKPM=B-1234\"')

	return parser.parse_args()


def main():
	args = parse_args()
	recs = SeqIO.parse(args.input, 'fasta')
	i = 1
	base_name = os.path.basename(args.input)
	if not args.locus_tag:
		locus_tag = os.path.splitext(base_name)[0]
	else:
		locus_tag = args.locus_tag
	out_name = (base_name, args.output)[bool(args.output)]
	
	if(os.path.exists(out_name)):
		print('Error! The file already exists!')
		exit(1)
	
	with open(out_name, 'w') as oh:
		for rec in recs:
			reclen = len(rec)
			if reclen < args.min_len:
				continue
			rec.id = f'{locus_tag}_{i}'
			rec.description = f'[{args.taxonomy } len={reclen}]'
			
			'''Unused regexp for coverage/depth
			sp = re.match(r'NODE_\d+_length_\d+_cov_([\d\.]+)', rec.id)
			uni = re.match(r'\d+ length=\d+ depth=([\d\.]+)x', rec.id)
			if sp:
				cov = sp.group(1)
			elif uni:
				depth = uni.group(1)
			else'''
			SeqIO.write(rec, oh, 'fasta')
			i += 1

if __name__ == '__main__':
	main()
