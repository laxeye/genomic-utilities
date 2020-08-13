#!/usr/bin/env python3
import dendropy
import sys
import io
import os
import subprocess
import argparse
import matplotlib.pyplot as plt
from matplotlib import cm
from numpy import array, arange
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list, fcluster, leaders
import logging


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
	"""Returns genome distances in form of dict"""
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


def ltmToDistDict(filename, is_aairb):
	"""Returns genome distances in form of dict"""
	try:
		lines = open(filename, 'r').readlines()
	except Exception as e:
		logger.error("Error: Impossible to read %s", filename)
		raise e
	lines.pop(0) # First line contains only number of taxa (rows)
	ids = []
	dists = {}

	for line in lines:
		if len(line) == 0:
			break
		l = line.rstrip("\r\n\t ").split('\t')
		name = l.pop(0).split("/")[-1]
		ids.append(name)
		dists[name] = {}
		for i, v in enumerate(l):
			dists[name][ids[i]] = aaiToDistance(v) if is_aairb else aniToDistance(v)

	return dists


def distDictToList(d):
	"""Create two-dimensional list from dict"""
	keys = list(d)
	ani_list = []

	for x in range(len(keys)):
		s = []
		for y in range(len(keys)):
			value = 0
			if x != y:
				(i, j) = (y, x) if x < y else (x, y)
				"""For compatibility with LT-matrices"""
				if ( keys[i] in d and keys[j] in d[keys[i]] ):
					value = d[keys[i]][keys[j]] 
				else:
					value = ANI_NA_DISTANCE
			s.append(value)
		ani_list.append(s)

	return ani_list


def printDistAsCSV(dists):
	"""Create in-memory file object"""
	inmem = io.StringIO()
	keys = list(dists)
	# 1st field [0,0] must be empty for Dendropy compatibility
	print("," + ",".join(keys), file=inmem)
	ani_list = distDictToList(dists)

	for i, key in enumerate(keys):
		t = list(map(str, ani_list[i]))
		t.insert(0, key)
		print(",".join(t), file=inmem)

	return inmem


def plotDendrogram(lm, names, prefix):
	logger.info("Plotting dendrogram")
	plt.figure(figsize=(10, 6), dpi=PNG_DPI)
	plt.title("Distance-derived dendrogram")
	dendrogram(lm, orientation='top', labels=names, distance_sort='descending', show_leaf_counts=True)
	plt.savefig("%s.dendrogram.png" % prefix, dpi=PNG_DPI)
	plt.savefig("%s.dendrogram.svg" % prefix)


def clusterGenomes(dist_dict, args):
	#Set colors for heatmap
	COLOR_MAP = cm.RdGy

	logger.info("Clustering genomes")
	dist_arr = array( distDictToList(dist_dict) )
	names = array( list(dist_dict) )
	lm = linkage(dist_arr, method="median", optimal_ordering=False)

	if args.plot_dendrogram:
		plotDendrogram(lm, names, args.prefix)

	if args.print_clusters:
		cluster_distance_t = 0.05
		cluster_ids = fcluster(lm, t=cluster_distance_t, criterion="distance")
		clustered_genomes = list(zip(names, cluster_ids))
		nodes, cluster_leader_ids = leaders(lm, cluster_ids)
		logger.debug("Leader node and cluster id:")
		logger.debug(list(zip(nodes, cluster_leader_ids)))
		with open("%s.repr.clstr" % args.prefix, 'w') as handle:
			for x in sorted(cluster_leader_ids):
				g, c = [y for y in clustered_genomes if y[1] == x][0]
				handle.write("%s\t%s\n" % (g, c))
		with open("%s.clstr" % args.prefix, 'w') as handle:
			for genome, cluster_id in clustered_genomes:
				handle.write("%s\t%s\n" % (genome, cluster_id))

	order = leaves_list(lm)
	dist_arr = dist_arr[order, ]
	dist_arr = dist_arr[:, order]
	names = names[order]

	"""Plot heatmap"""
	fig, ax = plt.subplots()
	ax.set_yticks( arange( len(dist_dict) ) )
	ax.set_xticks( arange( len(dist_dict) ) )
	ax.set_yticklabels(names, fontdict={'fontsize':'xx-small'})
	ax.set_xticklabels(names, fontdict={'fontsize':'xx-small'})
	plt.setp(ax.get_yticklabels(), rotation=0, ha="right", rotation_mode="anchor")
	plt.setp(ax.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor")

	im = ax.imshow(dist_arr, cmap=COLOR_MAP)
	ax.set_title("Distance between genomes")
	plt.colorbar(im)
	plt.savefig("%s.heatmap.png" % args.prefix, dpi=PNG_DPI)
	plt.savefig("%s.heatmap.svg" % args.prefix)


def pdmFromDistDict(dist_dict, prefix):
	#Returns phylogenetic distance matrix object from Dendropy
	outputtmp = printDistAsCSV(dist_dict)
	outputtmp.seek(0) # Return to the first line
	print(outputtmp.read(), file = open("%s.distances.csv" % prefix, 'w'))
	
	outputtmp.seek(0) # Return to the first line
	pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=outputtmp, delimiter=",")
	outputtmp.close()
	
	return pdm


def checkExec(binary):
	try:
		return subprocess.run([binary,'-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE).returncode is not None
	except Exception as e:
		raise e


def ani_from_report(filename):
	"""Return ANI value from dnadiff report"""
	with open(filename) as f:
		for line in f:
			l = line.rstrip().split()
			if len(l) > 0 and l[0] == "AvgIdentity":
				return (l[1])
		else:
			return "NA"


def get_dnadiff_report(path1, path2, prefix="tmp"):
	"""Run dnadiff and retrun path to report"""
	cmd = ["dnadiff", "-p", prefix, path1, path2]
	r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	if r.returncode == 0:
		return "%s.report" % prefix
	else:
		raise Exception("An error acquired during dnadiff execution!")


def listToFile(l, prefix):
	name = "%s.lst" % prefix
	with open(name,'w') as f:
		f.write("\n".join(l))
	return name


def parse_args():
	parser = argparse.ArgumentParser(description="Phylogenetic tree inference and heatmap drawing from genomic distances (based on ANI or AAI).")

	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument("-l", "--low-triangular-matrix", help="Low triangular matrix of distance values (ANI or AAI)")
	group.add_argument("-t", "--table", help="Tab separated table of distance values")
	group.add_argument("--anirb", help="Calculate ANI with ani.rb (slow). --input_list/--input_dir required", action="store_true")
	group.add_argument("--mummer",help="Calculate ANI with mummer. --input_list/--input_dir required", action="store_true")
	group.add_argument("--fastani", help="Calculate ANI with fastANI (fast). --input_list/--input_dir required", action="store_true")
	group.add_argument("--aairb", help="Calculate AAI with aai.rb (slow). --input_list/--input_dir required", action="store_true")

	ingroup = parser.add_mutually_exclusive_group()
	ingroup.add_argument("--input-list", help="List of genomes with full path for distance calculation")
	ingroup.add_argument("--input-dir", help="Path to directory containig genomes for distance calculation")
	parser.add_argument("-x", "--extension", help="Fasta files extension, e.g. fna (default), fa, fasta", default="fna")

	parser.add_argument("-p", "--prefix", help="Prefix for output files", default="newprefix")
	parser.add_argument("-m", "--tree-method", help="Phylogenetic tree inference method (default UPGMA)",
		default="UPGMA", choices=["UPGMA", "NJ", "both", "none"])
	parser.add_argument("-H", "--heatmap", help="Draw a heatmap", action="store_true")
	parser.add_argument("-A", "--ascii-tree", help="Draw ASCII tree to stdout", action="store_true")
	parser.add_argument("-d", "--plot-dendrogram", help="Plot a dendrogram", action="store_true")
	parser.add_argument("-c", "--print-clusters", help="Print genomic clusters", action="store_true")
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
		logger.critical("Unknown executable provided!")
		raise RuntimeError("Unknown executable provided!")


def calculateAI(args):
	binary = what_to_run(args)
	if not checkExec(binary):
		logger.critical("Failed to run \'%s\'", binary)
		sys.exit(1)

	file_results = "%s.tsv" % args.prefix

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
		logger.critical("Genomes are not provided. Exiting.")
		sys.exit(1)
	
	if len(file_list) < 2:
		logger.critical("Genome count too low: %s", len(file_list))
		logger.critical("Please check file extensions and path to files.")
		sys.exit(1)

	if args.anirb or args.aairb or args.mummer:
		results = []
		logger.info("Running %s", binary)
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
		with open(file_results, "w") as handle:
			handle.write("\n".join(results))
	
	if args.fastani:
		logger.info("Running fastANI")
		file_list_name = listToFile(file_list, args.prefix)
		cmd = ["fastANI", "--ql", file_list_name, "--rl", file_list_name, "-o", file_results, "-t", str(args.threads)]
		rc = subprocess.run(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE).returncode
		if rc != 0:
			logger.critical("fastANI didn't work correctly.")
			raise RuntimeError("fastANI error")

	return tableToDistDict(file_results, args.aairb)


def main():
	args = parse_args()
	global logger
	logger = logging.getLogger("main")
	logger.setLevel(logging.DEBUG)
	formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
	ch = logging.StreamHandler()
	ch.setLevel(logging.INFO)
	ch.setFormatter(formatter)
	logger.addHandler(ch)
	fh = logging.FileHandler("%s.log" % args.prefix)
	fh.setLevel(logging.DEBUG)
	fh.setFormatter(formatter)
	logger.addHandler(fh)
	logger.info("Analysis started")

	if (args.anirb or args.fastani or args.mummer or args.aairb):
		dist_dict = calculateAI(args)

	# Parse input files
	if args.low_triangular_matrix:
		dist_dict = ltmToDistDict(args.low_triangular_matrix, args.aairb)
	if args.table:
		dist_dict = tableToDistDict(args.table, args.aairb)

	# Clustering genomes, draw a heatmap and dendrogram with pyplot
	clusterGenomes(dist_dict, args)

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

	logger.info("Analysis finished")

if __name__ == '__main__':
	main()
