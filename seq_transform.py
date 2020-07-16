#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import os.path

def parse_args():
	parser = argparse.ArgumentParser(description="Reverse and/or complement.")

	parser.add_argument("-f", "--file", help="Input FASTA file")
	parser.add_argument("-s", "--sequence", help="Input sequence")
	parser.add_argument("-n", "--name", help="Sequence name", default="seq")
	parser.add_argument("-m", "--mode", default='all',
		help="Magic mode: all, rev, compl, revcompl",
		choices=['all', 'rev', 'compl', 'rebcompl'])
	parser.add_argument("-k", "--keep", action="store_true", help="Keep input sequence")

	args = parser.parse_args()

	if args.file and not os.path.exists(args.file):
		print(f"Error! File not found: '{args.file}'")
		raise FileNotFoundError(args.file)

	return args


def compl(s):
	if s == 'A':
		return 'T'
	elif s == 'C':
		return 'G'
	elif s == 'G':
		return 'C'
	elif s == 'T':
		return 'A'
	elif s == 'a':
		return 't'
	elif s == 'c':
		return 'g'
	elif s == 'g':
		return 'c'
	elif s == 't':
		return 'a'
	else:
		return 'N'


def rev_compl(seq):
	return ''.join(list(map(compl, seq[::-1])))


def rev_seq(seq):
	return seq[::-1]


def compl_seq(seq):
	return ''.join(list(map(compl, seq)))


def main():
	args = parse_args()
	seqlist = []

	if args.sequence:
		seqlist.append((args.name, args.sequence))

	if args.file:
		with open(os.path.abspath(args.file)) as fastafile:
			try:
				seqlist = list(SimpleFastaParser(fastafile))
			except Exception as e:
				raise e

	for (title, sequence) in seqlist:
		if args.keep:
			print(f">{title}\n{sequence}")
		if args.mode in ["all", "revcompl"]:
			print(f">{title} Reverse complement\n{rev_compl(sequence)}")
		if args.mode in ["all", "rev"]:
			print(f">{title} Reverse\n{rev_seq(sequence)}")
		if args.mode in ["all", "compl"]:
			print(f">{title} Complement\n{compl_seq(sequence)}")


if __name__ == '__main__':
	main()
