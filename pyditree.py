import dendropy
import sys
import io
import argparse
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
#from dendropy.datamodel import taxonmodel as taxonmodel
#from dendropy.calculate import phylogeneticdistance as phylodist
#from dendropy.utility import container

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

def heatmapFromDist(dist_dict):
	dist_arr = np.array( distDictToArray(dist_dict) )
	names = np.array( list(dist_dict) )

	fig, ax = plt.subplots()
	ax.set_yticks(np.arange(len(dist_dict)))

	if args.sort_heatmap:
		ordnung = np.argsort( np.sum(dist_arr, 0) )
		dist_arr = dist_arr[ordnung, ]
		dist_arr = dist_arr[:, ordnung ]
		names = names[ordnung]

	ax.set_yticklabels(names)
	plt.setp(ax.get_yticklabels(), rotation=30, ha="right", rotation_mode="anchor")

	im = ax.imshow(dist_arr, cmap=cm.RdGy)
	ax.set_title("ANI distance between genomes")
	plt.colorbar(im)
	plt.savefig("%s.png" % args.prefix, dpi=180)

def pdmFromDistDict(dist_dict):
	outputtmp = printDistAsCSV(dist_dict)
	
	outputtmp.seek(0) # Return to the first line
	print(outputtmp.read(), file = open("%s.dp.distances.csv" % args.prefix, 'w'))
	
	outputtmp.seek(0) # Return to the first line
	pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=outputtmp, delimiter=",")
	
	outputtmp.close()
	
	return pdm


parser = argparse.ArgumentParser()
parser.add_argument("-l", "--low_triangular_matrix", help="Low triangular matrix of ANI values")
parser.add_argument("-t", "--ani_table", help="Tab separated table of ANI values")
parser.add_argument("-p", "--prefix", help="Tab separated table of ANI values", default="newprefix")
parser.add_argument("-m", "--mode", help="Tree inference method: UPGMA (default), NJ or none", default="UPGMA")
parser.add_argument("-H", "--heatmap", help="Draw a heatmap", action="store_true")
parser.add_argument("-s", "--sort_heatmap", help="Sort heatmap rows and columns by cummulative distance", action="store_true")
parser.add_argument("-A", "--ascii_tree", help="Draw ASCII tree to stdout", action="store_true")
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

if args.mode == "UPGMA":
	upgma_tree = pdm.upgma_tree()
	print(upgma_tree.as_string(schema='newick', suppress_rooting = True), file = open("%s.upgma.unrooted.nwk" % args.prefix, 'w'))
	if args.ascii_tree:
		print(upgma_tree.as_ascii_plot(plot_metric='length'))
	upgma_tree.reroot_at_midpoint(suppress_unifurcations = False)
	print(upgma_tree.as_string(schema='newick', suppress_rooting = True), file = open("%s.upgma.rooted.nwk" % args.prefix, 'w'))
elif ars.mode == "NJ":
	nj_tree = pdm.nj_tree()
	print(nj_tree.as_string(schema='newick', suppress_rooting = True), file = open("%s.nj.unrooted.nwk" % args.prefix, 'w'))
	if args.ascii_tree:
		print(nj_tree.as_ascii_plot(plot_metric='length'))
	nj_tree.reroot_at_midpoint(suppress_unifurcations = False)
	print(nj_tree.as_string(schema='newick', suppress_rooting = True), file = open("%s.nj.rooted.nwk" % args.prefix, 'w'))
