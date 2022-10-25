#!/usr/bin/env python
from Bio import Entrez
import argparse
import re
import sqlite3
import sys


def assembly_info(assembly):
	result = dict()
	assembly_id = Entrez.read(Entrez.esearch(db='assembly', term=assembly))['IdList']
	if len(assembly_id) == 0:
		return None
	summary = Entrez.read(
		Entrez.esummary(db='assembly', id=assembly_id[0]), validate=False
	)['DocumentSummarySet']['DocumentSummary'][0]
	'''Summary example:
	{'DocumentSummarySet': DictElement({'DocumentSummary': [DictElement({'RsUid': '255868'
	'Taxid': '456320', 
	'Organism': 'Methanococcus voltae A3 (euryarchaeotes)',
	'SpeciesTaxid': '2188',
	'SpeciesName': 'Methanococcus voltae'
	'''

	if summary:
		result['id'] = assembly_id[0]
		result['SpeciesTaxid'] = summary['SpeciesTaxid']
		result['SpeciesName'] = summary['SpeciesName']
		result['Taxid'] = summary['Taxid']
		result['Organism'] = summary['Organism']
		return result
	else:
		print(f'Error with {assembly}')
		return None


def main(args):
	Entrez.email = 'email@domain.edu'
	conn = sqlite3.connect(args.database)
	c = conn.cursor()
	create_cmd = '''
		CREATE TABLE IF NOT EXISTS assemblies
		(acc text, id text, taxid text, organism text,
		speciestaxid text, species text)
	'''
	conn.execute(create_cmd)

	insert_cmd = '''
		INSERT INTO assemblies (acc, id, taxid, organism, speciestaxid, species)
		VALUES (?, ?, ?, ?, ?, ?)
	'''
	
	rgx = re.compile('GC[AF]_\d+\.\d')

	try:
		h = open(args.input, 'r')
	except Exception as e:
		raise e

	for x in h.readlines():
		x = x.rstrip()
		m = rgx.search(x)
		if m:
			acc = m.group()
		else:
			continue

		# Check if accession is already in the DB
		c.execute('SELECT EXISTS(SELECT 1 FROM assemblies where acc=?)', (acc,))
		if c.fetchone()[0] == 1:
			if args.print:
				res = c.execute('SELECT * FROM assemblies where acc=?', (acc,))
				print('\t'.join(res.fetchone()))
		else:
			r = assembly_info(acc)
			if r:
				if args.print:
					print('\t'.join([acc] + list(r.values())))
				c.execute(insert_cmd, (
					acc, r['id'], r['Taxid'], r['Organism'],
					r['SpeciesTaxid'], r['SpeciesName']
				))
				conn.commit()

	h.close()


def parse_args():
	parser = argparse.ArgumentParser(
		description='Find and store taxonomy of NCBI Assembly accessions.'
		)
	parser.add_argument('-d', '--database', required=True,
		help='SQLite3 database to store the data.')
	parser.add_argument('-p', '--print', action='store_true',
		help='Print found taxonomies to stdout.')
	
	in_group = parser.add_mutually_exclusive_group(required=True)
	in_group.add_argument('-i', '--input',
		help='List of accessions, one per line.')
	in_group.add_argument('-q', '--query',
		help='Single accession.')
	
	return parser.parse_args()


if __name__ == '__main__':
	main(parse_args())
