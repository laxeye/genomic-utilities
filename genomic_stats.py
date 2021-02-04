#!/usr/bin/env python3
'''Calculate statistics of genome assembly including:
GC content, sequence count, total length, maximum sequence length,
N50, L50, N90, L90 and user-provided Nx and Lx.
'''
import sys
import json
import argparse
import sqlite3
import os
import gzip
import bz2
from Bio import SeqIO
from Bio.SeqUtils import GC


def get_metric(lengths, value=50):
	'''Returns a tuple containing NX and LX metric'''
	l_total = sum(lengths)
	metric = 0.01 * value
	l_sum = 0
	for i, x in enumerate(lengths, 1):
		l_sum += x
		if l_sum >= l_total * metric:
			return x, i
	return None, None


def assembly_stats(genome, extra=None):
	'''Returns a dict containing genome stats.

	Additional metric extra (int) may be provided to calculate Nx and Lx.
	GC content is rounded to two decimal places
	'''
	stats = dict()
	extension = genome.split(".")[-1]
	if extension == 'gz':
		try:
			handle = gzip.open(genome, 'rt')
		except Exception as e:
			raise e
	elif extension == 'bz2':
		try:
			handle = bz2.open(genome, 'rt')
		except Exception as e:
			raise e
	else:
		try:
			handle = open(genome, 'r')
		except Exception as e:
			raise e

	seq_records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	lengths = sorted([len(record.seq) for record in seq_records], reverse=True)
	stats['Filename'] = os.path.basename(genome)
	stats['GC'] = "{:.2f}".format(round(GC(''.join(map(lambda it: str(it.seq), seq_records))), 2))
	stats['Sequence count'] = len(lengths)
	stats['Total length'] = sum(lengths)
	stats['Max length'] = lengths[0]
	stats['N50'], stats['L50'] = get_metric(lengths, 50)
	stats['N90'], stats['L90'] = get_metric(lengths, 90)
	if extra:
		stats[f'N{extra}'], stats[f'L{extra}'] = get_metric(lengths, extra)
	return stats


def write_assembly_stats(stats, args):
	'''Writes assembly stats to file or stdout'''
	prefix = args.output_prefix
	s_format = args.format if args.format else "human"

	if prefix:
		ext_dict = {"human": "txt", "json": "json", "table": "tsv"}
		filename = f"{prefix}.{ext_dict[s_format]}"
		mode = 'w' if args.overwrite else 'a'
		try:
			dest = open(filename, mode)
		except Exception as e:
			raise e
	else:
		dest = sys.stdout

	if s_format == "human":
		for genome in stats:
			for k, v in genome.items():
				print(f"{k}\t{v}", file=dest)
	if s_format == "json":
		print(json.dumps(stats, indent=4), file=dest)
	if s_format == "table":
		header = "\t".join(stats[0].keys())
		print(f"#{header}", file=dest)
		for genome in stats:
			data = "\t".join(list(map(str, genome.values())))
			print(f"{data}", file=dest)

	if prefix:
		dest.close()


def write_to_db(stats, db):
	'''Write assembly stats to SQLite database'''
	con = sqlite3.connect(db)
	cur = con.cursor()
	cur.execute("CREATE TABLE IF NOT EXISTS stats (filename, GC, seq_count, total_length, max_length, N50, L50, N90, L90);")
	for genome in stats:
		cur.execute("INSERT INTO stats values (?,?,?,?,?,?,?,?,?)", tuple(genome.values()))
	con.commit()
	con.close()


def parse_args():
	'''Retruns parsed args'''
	parser = argparse.ArgumentParser(description="Genome statistics")
	parser.add_argument("-i", "--input", required=True,
		help="Input file in FASTA format or directory.")
	parser.add_argument("-o", "--output-prefix",
		help="Output prefix. Skip it to print to stdout. "
		+ "If file already exists output will be appended.")
	parser.add_argument("-w", "--overwrite", action="store_true",
		help="Overwrite output file.")
	parser.add_argument("-f", "--format", default="human",
		choices=["human", "json", "table"], help="Output format: "
		+ "human (human-fiendly, default), json, table (tab-delimited).")
	parser.add_argument("-X", "--additional-metric", type=int,
		help="Additional metric to calculate Nx and Lx. "
		+ "Integer between 1 and 99")
	parser.add_argument("-s", "--sqlite-db",
		help="Path to sqlite database to store data.")
	parser.add_argument("-e", "--extension",
		help="File extension if directory was provided. "
		+ "By default any file will be processed.")
	parser.add_argument("-b", "--basename-as-prefix", action="store_true",
		help="Ouput file prefix will be taken from input file or directory.")

	args = parser.parse_args()

	args.input = os.path.abspath(args.input)
	if not os.path.exists(args.input):
		raise FileNotFoundError(args.input)
	if args.basename_as_prefix and not args.prefix:
		args.prefix = os.path.basename(args.input)
	if args.additional_metric and (
			args.additional_metric >= 100 or args.additional_metric <= 0):
		print(f"Illegal metric value: {args.additional_metric}.",
			file=sys.stderr)
		raise ValueError("Illegal additional metric value!")

	return args


def main():
	args = parse_args()
	file_list = []
	if os.path.isdir(args.input):
		for file in os.scandir(os.path.realpath(args.input)):
			if file.is_file():
				extension = file.name.split(".")[-1]
				if extension in ['gz', 'bz2']:
					extension = file.name.split(".")[-2]
				if (not args.extension
					or extension == args.extension
				):
					file_list.append(file.path)
	else:
		file_list.append(args.input)

	stats = [assembly_stats(genome, args.additional_metric) for genome in file_list]
	write_assembly_stats(stats, args)
	if args.sqlite_db:
		write_to_db(stats, args.sqlite_db)


if __name__ == '__main__':
	main()
