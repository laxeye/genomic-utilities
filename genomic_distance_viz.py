#!/usr/bin/env python3
import time
import dendropy
import sys
import io
import os
import subprocess
import argparse
import matplotlib.pyplot as plt
from matplotlib import cm
from numpy import array, arange
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list


#Set distance for unknown (low) ANI and AAI
ANI_NA_DISTANCE = 0.3
AAI_NA_DISTANCE = 0.8
#PNG resolution
PNG_DPI = 300


def aniToDistance(ani):
	return ANI_NA_DISTANCE if (ani == "NA") else ( 100 - float(ani) ) / 100.0


def aaiToDistance(ani):
	return AAI_NA_DISTANCE if (ani == "NA") else ( 100 - float(ani) ) / 100.0


def tableToDistDict(filename, is_aairb):
	dists = {}
	ids = []
	
	with open(filename, 'r') as infile:
		for line in infile:
			(query, ref, value) = line.rstrip('\n\r\t ').split("\t")[0:3]
			query = query.split("/")[-1]
			ref = ref.split("/")[-1]
			if query not in ids:
				dists[query] = {}
				ids.append(query)
			if ref not in ids:
				dists[ref] = {}
				ids.append(ref)

			value = aaiToDistance(value) if is_aairb else aniToDistance(value)
			dists[query][ref] = value
			dists[ref][query] = value

	return dists


def ltmToDistDict(filename):
	lines = open(filename, 'r').read().split("\n")
	lines.pop(0) # First line contains only number of taxa (rows)
	ids = []
	dists = {}

	for line in lines:
		if len(line) == 0: break
		l = line.rstrip("\r\n\t ").split('\t')
		name = l.pop(0).split("/")[-1]
		ids.append(name)
		dists[name] = {}
		for i in range(len(l)):
			#Supposing ANI
			dists[name][ids[i]] = aniToDistance(l[i])

	return dists


def distDictToList(d):
	#Creates two-dimensional list from dict
	keys = list(d)
	ani_list = []

	for x in range(len(keys)):
		s = []
		for y in range(len(keys)):
			value = 0
			if x != y:
				(i, j) = (y, x) if x < y else (x, y) # Compatibility with LT-matrices
				if ( keys[i] in d and keys[j] in d[keys[i]] ):
					value = d[keys[i]][keys[j]] 
				else:
					value = ANI_NA_DISTANCE
			s.append(value)
		ani_list.append(s)

	return ani_list


def printDistAsCSV(dists):
	inmem = io.StringIO()
	keys = list(dists)
	# 1st field [0,0] must be empty for Dendropy compatibility
	print("," + ",".join(keys), file = inmem)
	ani_list = distDictToList(dists)

	for x in range(len(keys)):
		t = list(map(lambda x: str(x), ani_list[x]))
		t.insert(0, keys[x])
		print(",".join(t), file = inmem)

	return inmem


def plotDendrogram(l, names, prefix):
	plt.figure(figsize=(10, 6), dpi=PNG_DPI)
	plt.title("Distance-derived dendrogram")
	dendrogram(l, orientation='top', labels=names, distance_sort='descending', show_leaf_counts=True)
	plt.savefig("%s.dendrogram.png" % prefix, dpi=PNG_DPI)
	plt.savefig("%s.dendrogram.svg" % prefix)


def heatmapFromDist(dist_dict, prefix, plot_dendrogram):
	#Set colors for heatmap
	COLOR_MAP = cm.RdGy

	dist_arr = array( distDictToList(dist_dict) )
	names = array( list(dist_dict) )
	l = linkage(dist_arr)

	if plot_dendrogram:
		plotDendrogram(l, names, prefix)

	fig, ax = plt.subplots()
	ax.set_yticks( arange( len(dist_dict) ) )
	ax.set_xticks( arange( len(dist_dict) ) )

	order = leaves_list(l)
	dist_arr = dist_arr[order, ]
	dist_arr = dist_arr[:, order]
	names = names[order]

	ax.set_yticklabels(names, fontdict={'fontsize':'xx-small'})
	ax.set_xticklabels(names, fontdict={'fontsize':'xx-small'})
	plt.setp(ax.get_yticklabels(), rotation=0, ha="right", rotation_mode="anchor")
	plt.setp(ax.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor")

	im = ax.imshow(dist_arr, cmap=COLOR_MAP)
	ax.set_title("Distance between genomes")
	plt.colorbar(im)
	plt.savefig("%s.heatmap.png" % prefix, dpi=PNG_DPI)
	plt.savefig("%s.heatmap.svg" % prefix)


def pdmFromDistDict(dist_dict, prefix):
	#Returns phylogenetic distance matrix object from Dendropy
	outputtmp = printDistAsCSV(dist_dict)
	outputtmp.seek(0) # Return to the first line
	print(outputtmp.read(), file = open("%s.distances.csv" % prefix, 'w'))
	
	outputtmp.seek(0) # Return to the first line
	pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=outputtmp, delimiter=",")
	outputtmp.close()
	
	return pdm


def checkExec(bin):
	try:
		return subprocess.run([bin,'-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).returncode is not None
	except Exception as e:
		raise e

#dnadiff .report parsing
def ani_from_report(filename):
	with open(filename) as f:
		for line in f:
			l = line.rstrip().split()
			if len(l) > 0 and l[0] == "AvgIdentity":
				return (l[1])

#dnadiff run
def get_dnadiff_report(path1, path2, prefix="tmp"):
	cmd = ["dnadiff", "-p", prefix, path1, path2]
	r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if r.returncode == 0:
		return f"{prefix}.report"
	else:
		raise Exception("An error acquired during dnadiff execution!")


def listToFile(l, prefix):
	name = "%s.lst" % prefix
	with open(name,'w') as f:
		f.write("\n".join(l))
	return name


def parse_args():
	parser = argparse.ArgumentParser(description = "Phylogenetic tree inference and heatmap drawing from genomic distances (based on ANI or AAI).")

	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument("-l", "--low_triangular_matrix", help="Low triangular matrix of distance values (ANI or AAI)")
	group.add_argument("-t", "--table", help="Tab separated table of distance values")
	group.add_argument("--anirb", help="Calculate ANI with ani.rb (slow). --input_list/--input_dir required", action="store_true")
	group.add_argument("--mummer",help="Calculate ANI with mummer. --input_list/--input_dir required", action="store_true")
	group.add_argument("--fastani", help="Calculate ANI with fastANI (fast). --input_list/--input_dir required", action="store_true")
	group.add_argument("--aairb", help="Calculate AAI with aai.rb (slow). --input_list/--input_dir required", action="store_true")

	ingroup = parser.add_mutually_exclusive_group()
	ingroup.add_argument("--input_list", help="List of genomes with full path for distance calculation")
	ingroup.add_argument("--input_dir", help="Path to directory containig genomes for distance calculation")
	parser.add_argument("-x", "--extension", help="Fasta files extension, e.g. fna (default), fa, fasta", default="fna")

	parser.add_argument("-p", "--prefix", help="Prefix for output files", default="newprefix")
	parser.add_argument("-m", "--tree_method", help="Phylogenetic tree inference method (default UPGMA)", default="UPGMA", choices=["UPGMA", "NJ", "both", "none"])
	parser.add_argument("-H", "--heatmap", help="Draw a heatmap", action="store_true")
	parser.add_argument("-A", "--ascii_tree", help="Draw ASCII tree to stdout", action="store_true")
	parser.add_argument("-d", "--plot_dendrogram", help="Plot a dendrogram", action="store_true")
	parser.add_argument("--reroot", help="Reroot tree at midpoint. May cause errors or incorrect trees", action="store_true")
	parser.add_argument("--threads", help="Number of CPU threads (where possible)", type=int, default=1)

	return parser.parse_args()


def what_to_run(args):
	if args.anirb:
		return "ani.rb"
	elif args.aairb:
		return "aai.rb"
	elif args.fastani:
		return "fastANI"
	elif args.mummer:
		return "dnadiff"
	else:
		raise Exception("Unknown executable provided!")
		return None

def main():
	args = parse_args()

	if (args.anirb or args.fastani or args.mummer or args.aairb):
		binary = what_to_run(args)
		if checkExec(binary) == False:
			print("Failed to run \'%s\'" % binary)
			exit(1)
		file_list =[]
		if args.input_list:
			f = open(args.input_list)
			for l in f.readlines():
				l = l.strip()
				if (os.path.exists(l) and l.split(".")[-1] == args.extension):
					file_list.append(os.path.realpath(l))
		elif args.input_dir:
			for file in os.scandir(os.path.realpath(args.input_dir)):
				if (file.is_file() and file.name.split(".")[-1] == args.extension):
					file_list.append(file.path)
		else:
			print("Genomes are not provided. Exiting.")
			exit(1)
		
		if len(file_list) < 2:
			print(f"Genome count too low: {len(file_list)}.\n" + 
				"Please check file extensions and path to files.")
			exit(1)

		if args.anirb or args.aairb or args.mummer:
			results = []
			print(f"Running {binary}", file = sys.stderr)
			#Generate list of pairs (Genome1, Genome2) for A*I calculating
			pairs = [(file_list[i], file_list[j]) for i in range(len(file_list)) for j in range(i,len(file_list)) if i != j]
			for (i, j) in pairs:
				if args.mummer:
					report = get_dnadiff_report(i,j)
					r = ani_from_report(report)
				else:
					cmd = [binary, "-1", i, "-2", j, "-a", "-q", "-t", str(args.threads)]
					if args.aairb:
						cmd = cmd + ['-p', 'diamond']
					r = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.DEVNULL, text = True).stdout.strip()
				if len(r) > 0:
					results.append("\t".join([os.path.basename(i), os.path.basename(j), r]))
			for file in file_list:
				results.append("\t".join([os.path.basename(file), os.path.basename(file), str(100)]))
			file_results = "%s.tsv" % args.prefix
			print("\n".join(results), file = open(file_results, "w"))
		
		if args.fastani:
			print("Running fastANI", file = sys.stderr)
			file_list_name = listToFile(file_list, args.prefix)
			file_results = "%s.tsv" % args.prefix
			cmd = ["fastANI", "--ql", file_list_name, "--rl", file_list_name, "-o", file_results, "-t", str(args.threads)]
			rc = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE).returncode
			if rc != 0:
				print("fastANI didn't work correctly.")
				exit(1)

		dist_dict = tableToDistDict(file_results, args.aairb)

	# Parse input files
	if args.low_triangular_matrix:
		dist_dict = ltmToDistDict(args.low_triangular_matrix)
	if args.table:
		dist_dict = tableToDistDict(args.table, args.aairb)

	# Draw a heatmap with pyplot
	if args.heatmap:
		heatmapFromDist(dist_dict, args.prefix, args.plot_dendrogram)

	# Calculate Dendropy compatible phylogenetic distance matrix
	pdm = pdmFromDistDict(dist_dict, args.prefix)

	if (args.tree_method == "UPGMA" or args.tree_method == "both"):
		upgma_tree = pdm.upgma_tree()
		print(upgma_tree.as_string(schema='newick', suppress_rooting = True), file = open("%s.upgma.unrooted.nwk" % args.prefix, 'w'))
		if args.ascii_tree:
			print(upgma_tree.as_ascii_plot(plot_metric='length'), file = open("%s.upgma.unrooted.txt" % args.prefix, 'w'))
		if args.reroot:
			upgma_tree.reroot_at_midpoint(suppress_unifurcations = False)
			print(upgma_tree.as_string(schema='newick', suppress_rooting = True), file = open("%s.upgma.rooted.nwk" % args.prefix, 'w'))
			if args.ascii_tree:
				print(upgma_tree.as_ascii_plot(plot_metric='length'), file = open("%s.upgma.rooted.txt" % args.prefix, 'w'))
		
	if (args.tree_method == "NJ" or args.tree_method == "both"):
		nj_tree = pdm.nj_tree()
		print(nj_tree.as_string(schema='newick', suppress_rooting = True), file = open("%s.nj.unrooted.nwk" % args.prefix, 'w'))
		if args.ascii_tree:
			print(nj_tree.as_ascii_plot(plot_metric='length'), file = open("%s.nj.unrooted.txt" % args.prefix, 'w'))
		if args.reroot:
			nj_tree.reroot_at_midpoint(suppress_unifurcations = False)
			print(nj_tree.as_string(schema='newick', suppress_rooting = True), file = open("%s.nj.rooted.nwk" % args.prefix, 'w'))
			if args.ascii_tree:
				print(nj_tree.as_ascii_plot(plot_metric='length'), file = open("%s.nj.rooted.txt" % args.prefix, 'w'))

if __name__ == '__main__':
	main()
