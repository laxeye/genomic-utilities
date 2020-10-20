#!/usr/bin/env python3
'''Calculates statistics of genome assembly including:
Sequence count, total length, maximum sequence length,
N50, L50, N90, L90 and user-provided Nx and Lx.
'''
import sys
import json
import argparse
from Bio import SeqIO
import sqlite3
import os


def get_N_L_metric(lengths, value=50):
	'''Returns a tuple containing NX and LX metric'''
	l_total = sum(lengths)
	metric = 0.01 * value
	l_sum = 0
	for i, x in enumerate(lengths, 1):
		l_sum += x
		if l_sum >= l_total * metric:
			return x, i
	return None, None


def assembly_stats(genome, x=None):
	'''Returns a dict containing genome stats.

	Additional metric x (int) may be provided to calculate Nx and Lx.
	'''
	stats = dict()
	seq_records = SeqIO.parse(genome, "fasta")
	lengths = sorted([len(x.seq) for x in seq_records], reverse=True)
	stats['Sequence count'] = len(lengths)
	stats['Total length'] = sum(lengths)
	stats['Max length'] = lengths[0]
	stats['N50'], stats['L50'] = get_N_L_metric(lengths, 50)
	stats['N90'], stats['L90'] = get_N_L_metric(lengths, 90)
	if x:
		stats[f'N{x}'], stats[f'L{x}'] = get_N_L_metric(lengths, x)
	return stats


def write_assembly_stats(stats, prefix=None, s_format="human"):
	'''Writes assembly stats to file or stdout'''
	if prefix:
		ext_dict = {"human": "txt", "json": "json", "table": "tsv"}
		filename = f"{prefix}.{ext_dict[s_format]}"
		try:
			dest = open(filename, 'w')
		except Exception as e:
			raise e
	else:
		dest = sys.stdout

	if s_format == "human":
		for k, v in stats.items():
			print(f"{k}\t{v}", file=dest)
	if s_format == "json":
		print(json.dumps(stats, indent=4), file=dest)
	if s_format == "table":
		header = "\t".join(stats.keys())
		data = "\t".join(list(map(str, stats.values())))
		print(f"#{header}", file=dest)
		print(f"{data}", file=dest)


def write_to_db(stats, db, filename):
	'''Write assembly stats to SQLite database'''
	con = sqlite3.connect(db)
	cur = con.cursor()
	filename = os.path.basename(filename)
	cur.execute("CREATE TABLE IF NOT EXISTS stats (genome, seq_count, total_length, max_length, N50, L50, N90, L90);")
	cur.execute("INSERT INTO stats values (?,?,?,?,?,?,?,?)", (filename ,*stats.values()))
	con.commit()
	con.close()



def parse_args():
	'''Retruns parsed args'''
	parser = argparse.ArgumentParser(description="Genome statistics")
	parser.add_argument("-i", "--input", required=True,
		help="Input file in FASTA format.")
	parser.add_argument("-o", "--output-prefix",
		help="Output prefix. Skip it to print to stdout.")
	parser.add_argument("-f", "--format", default="human",
		choices=["human", "json", "table"], help="Output format: "
		+ "human (human-fiendly, default), json, table (tab-delimited).")
	parser.add_argument("-X", "--additional-metric", type=int,
		help="Additional metric to calculate Nx and Lx. "
		+ "Integer between 1 and 99")
	parser.add_argument("-s", "--sqlite-db",
		help="Path to sqlite database to store data.")

	args = parser.parse_args()

	if args.additional_metric and (
			args.additional_metric >= 100 or args.additional_metric <= 0):
		print(f"Illegal metric value: {args.additional_metric}.",
			file=sys.stderr)
		raise ValueError("Illegal additional metric value!")

	return args


if __name__ == '__main__':
	args = parse_args()
	stats = assembly_stats(args.input, args.additional_metric)
	if args.sqlite_db:
		write_to_db(stats, args.sqlite_db, args.input)
	write_assembly_stats(
		stats, prefix=args.output_prefix, s_format=args.format)
