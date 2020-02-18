#!/usr/bin/env python3
import dendropy
import sys
import io
import argparse
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy.array, numpy.arange
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list

#Set distance for unknown (low) ANI
NA_DISTANCE = 0.3

def aniToDistance(ani):
	return NA_DISTANCE if (ani == "NA") else ( 100 - float(ani) ) / 100.0

def tableToDistDict(filename):
	infile = open(filename, 'r')
	dists = {}
	ids = []
	line = infile.readline()

	while line:
		line.rstrip("\n\r\t ")
		(query, ref, ani) = line.split("\t")[0:3]
		if query not in ids:
			dists[query] = {}
			ids.append(query)
		dists[query][ref] = aniToDistance(ani)
		line = infile.readline()

	return dists

def ltmToDistDict(filename):
	lines = open(filename, 'r').read().split("\n")
	lines.pop(0) # First line contains only number of taxa (rows)
	ids = []
	dists = {}

	for line in lines:
		if len(line) == 0: break
		l = line.rstrip("\r\n\t ").split('\t')
		name = l.pop(0)
		ids.append(name)
		dists[name] = {}
		for i in range(len(l)):
			dists[name][ids[i]] = aniToDistance(l[i])

	return dists

def distDictToArray(d):
	keys = list(d)
	a = []

	for x in range(len(keys)):
		s = []
		for y in range(len(keys)):
			v = 0
			if x != y:
				(i, j) = (y, x) if x < y else (x, y) # Compatibility with LT-matrices
				if ( keys[i] in d and keys[j] in d[keys[i]] ):
					v = d[keys[i]][keys[j]] 
				else:
					v = NA_DISTANCE # N/A tmp value
			s.append(v)
		a.append(s)

	return a

def printDistAsCSV(d):
	inmem = io.StringIO()
	keys = list(d)
	# 1st field [0,0] must be empty for Dendropy compatibility
	print("," + ",".join(keys), file = inmem)
	a = distDictToArray(d)

	for x in range(len(keys)):
		t = list(map(lambda x: str(x), a[x]))
		t.insert(0, keys[x])
		print(",".join(t), file = inmem)

	return inmem

def plotDendrogram(l, names):
	plt.figure(figsize=(10, 6), dpi=300)
	plt.title("ANI-derived dendrogram")
	dendrogram(l, orientation='top', labels=names, distance_sort='descending', show_leaf_counts=True)
	plt.savefig("%s.dendrogram.png" % args.prefix, dpi=300)
	plt.savefig("%s.dendrogram.svg" % args.prefix)

def heatmapFromDist(dist_dict):
	dist_arr = numpy.array( distDictToArray(dist_dict) )
	names = numpy.array( list(dist_dict) )
	l = linkage(dist_arr)
	fig, ax = plt.subplots()
	ax.set_yticks(numpy.arange(len(dist_dict)))

	order = leaves_list(l)
	dist_arr = dist_arr[order, ]
	dist_arr = dist_arr[:, order]
	names = names[order]

	ax.set_yticklabels(names,fontdict={'fontsize':'xx-small'})
	plt.setp(ax.get_yticklabels(), rotation=0, ha="right", rotation_mode="anchor")

	im = ax.imshow(dist_arr, cmap=cm.RdGy)
	ax.set_title("ANI distance between genomes")
	plt.colorbar(im)
	plt.savefig("%s.heatmap.png" % args.prefix, dpi=300)
	plt.savefig("%s.heatmap.svg" % args.prefix)

	if args.plot_dendrogram:
		plotDendrogram(l, names)

def pdmFromDistDict(dist_dict):
	outputtmp = printDistAsCSV(dist_dict)
	
	outputtmp.seek(0) # Return to the first line
	print(outputtmp.read(), file = open("%s.dp.distances.csv" % args.prefix, 'w'))
	
	outputtmp.seek(0) # Return to the first line
	pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=outputtmp, delimiter=",")
	
	outputtmp.close()
	
	return pdm


parser = argparse.ArgumentParser(description = "Phylogenetic tree inference and heatmap drawing from ANI-derived genomic distances.")
parser.add_argument("-l", "--low_triangular_matrix", help="Low triangular matrix of ANI values")
parser.add_argument("-t", "--ani_table", help="Tab separated table of ANI values")
parser.add_argument("-p", "--prefix", help="Tab separated table of ANI values", default="newprefix")
parser.add_argument("-m", "--mode", help="Tree inference method: UPGMA (default), NJ, both or none", default="UPGMA")
parser.add_argument("-H", "--heatmap", help="Draw a heatmap", action="store_true")
parser.add_argument("-A", "--ascii_tree", help="Draw ASCII tree to stdout", action="store_true")
parser.add_argument("-d", "--plot_dendrogram", help="Plot a dendrogram", action="store_true")
args = parser.parse_args()

# Parse input
if args.low_triangular_matrix:
	dist_dict = ltmToDistDict(args.low_triangular_matrix)
elif args.ani_table:
	dist_dict = tableToDistDict(args.ani_table)
else:
	print("You should provide ANI values.\nPlease run %s -h" % sys.argv[0])
	exit(1)

# Draw a heatmap with pyplot
if args.heatmap:
	heatmapFromDist(dist_dict)

# Calculate Dendropy compatible phylogenetic distance matrix
pdm = pdmFromDistDict(dist_dict)

if (args.mode == "UPGMA" or args.mode == "both"):
	upgma_tree = pdm.upgma_tree()
	print(upgma_tree.as_string(schema='newick', suppress_rooting = True), file = open("%s.upgma.unrooted.nwk" % args.prefix, 'w'))
	if args.ascii_tree:
		print(upgma_tree.as_ascii_plot(plot_metric='length'), file = open("%s.upgma.unrooted.txt" % args.prefix, 'w'))
	upgma_tree.reroot_at_midpoint(suppress_unifurcations = False)
	print(upgma_tree.as_string(schema='newick', suppress_rooting = True), file = open("%s.upgma.rooted.nwk" % args.prefix, 'w'))
if (args.mode == "NJ" or args.mode == "both"):
	nj_tree = pdm.nj_tree()
	print(nj_tree.as_string(schema='newick', suppress_rooting = True), file = open("%s.nj.unrooted.nwk" % args.prefix, 'w'))
	if args.ascii_tree:
		print(nj_tree.as_ascii_plot(plot_metric='length'), file = open("%s.nj.unrooted.txt" % args.prefix, 'w'))
	nj_tree.reroot_at_midpoint(suppress_unifurcations = False)
	print(nj_tree.as_string(schema='newick', suppress_rooting = True), file = open("%s.nj.rooted.nwk" % args.prefix, 'w'))
