# `genomic_distance_viz.py`

Phylogenetic tree inference and heatmap drawing from ANI (average nucleotide identity) or AAI (average aminoacid identity)-derived genomic distances.


## Usage

`genomic_distance_viz.py [-h] (-l LOW_TRIANGULAR_MATRIX | -t ANI_TABLE | --anirb | --aairb | --mummer | --fastani) [--input-list INPUT_LIST | --input-dir INPUT_DIR] [-d] [-x EXTENSION] [-p PREFIX] [-m {UPGMA,NJ,both,none}] [--threads THREADS] [-H] [-A] [-D] [--cluster-threshold CLUSTER_THRESHOLD] [--reroot]`


### Command-line options

```
-h, --help
Show this help message and exit

-l FILE, --low-triangular-matrix FILE
Low triangular matrix of ANI values

-t FILE, --table FILE
Tab separated table of similarity (ANI or AAI) between genomes

--anirb
Calculate ANI with ani.rb (should be installed separately) (slow). --input-list/--input-dir required

--fastani
Calculate ANI with fastANI (should be installed separately) (fast). --input_list/--input_dir required

--mummer
Calculate ANI with mummer. --input_list/--input_dir required

--input-list FILE
File with a list of full paths of genomes for ANI calculation

--input-dir DIRECTORY
Path (may be relative) to directory containig genomes for ANI calculation

-x EXTENSION, --extension EXTENSION
Fasta files extension, e.g. fna (default), fa, fasta

-d, --use-diamond
Use diamond instead of BLAST in aai.rb

-p PREFIX, --prefix PREFIX
Prefix for output files

-m {UPGMA,NJ,both,none}, --tree-method {UPGMA,NJ,both,none}
Phylogenetic tree inference method (default UPGMA)

-H, --heatmap
Draw a heatmap

-A, --ascii-tree
Draw ASCII tree to stdout

-D, --plot-dendrogram
Plot a dendrogram

-c, --print-clusters
Print genomic clusters (experimental)

--cluster-threshold CLUSTER_THRESHOLD
Threshold for genomic clusters output

--reroot
Reroot tree at midpoint. May cause errors or incorrect trees

--threads THREADS
Number of CPU threads (where possible)
```


### Input files

**Low triangular matrix** - matrix produced by [fastANI](https://github.com/ParBLiSS/FastANI) with "--matrix" option or any other software.

**ANI table** - tab-separated file of such structure: 

`Genome1    Genome2    Identity[   ...]`

May be produced by e.g. [ani.rb](https://github.com/lmrodriguezr/enveomics):

`for i in *fna; do for j in *fna; do echo -ne "$i\t$j\t"; ani.rb -q -a -1 $i -2 $j 2>/dev/null; done; done >> ani.rb.tsv`


Or by [fastANI](https://github.com/ParBLiSS/FastANI).


If you have **ani.rb**, **aai.rb**, **mummer4** or **fastANI** installed in Your environment you may use corresponding key to calculate genome identity for list of genomes using *--input_list* or folder containing genomes *--input_dir*. Genomes should have FASTA format, You may provide an extension using *--extension* option.


## Requirements

* dendropy
* matplotlib
* numpy
* scipy

May be installed with pip:

`pip install dendropy matplotlib numpy scipy `

or conda

`conda install dendropy matplotlib numpy scipy`

Optional:
* [fastANI](https://github.com/ParBLiSS/FastANI)
* [ani.rb](https://github.com/lmrodriguezr/enveomics)
* [aai.rb](https://github.com/lmrodriguezr/enveomics)
* [Mummer4](https://github.com/mummer4/mummer)

## Known bugs

AssertionError during midpoint rooting

`File "/somewhere/site-packages/dendropy/datamodel/treemodel.py", line 5076, in reroot_at_midpoint
assert break_on_node is not None or target_edge is not None
AssertionError`


# `genomic_stats.py`

basic genome statistics retreived with BioPython.

## Usage

`genomic_stats.py [-h] -i INPUT [-o OUTPUT_PREFIX] [-f {human,json,table}] [-X ADDITIONAL_METRIC]`

Arguments:

-i INPUT, --input INPUT  
Input file in FASTA format.

-o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX  
Output prefix. Skip it to print to stdout.

-f {human,json,table}, --format {human,json,table}  
Output format: human (human-fiendly, default), json,
table (tab-delimited).

-X ADDITIONAL_METRIC, --additional-metric ADDITIONAL_METRIC  
Additional metric to calculate Nx and Lx. Integer between 1 and 99

## Requirements

* Biopython

May be installed with pip:

`pip install biopython`

or conda

`conda install biopython`
