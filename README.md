# `ani_distance_viz.py`

Phylogenetic tree inference and heatmap drawing from ANI (average nucleotide identity)-derived genomic distances.


## Usage

`ani_distance_viz.py [-h] (-l LOW_TRIANGULAR_MATRIX | -t ANI_TABLE | --anirb | --fastani) [--input_list INPUT_LIST | --input_dir INPUT_DIR] [--extension EXTENSION] [-p PREFIX] [-m MODE] [-H] [-A] [-d] [--reroot]`


### Command-line options

-h, --help                              Show this help message and exit

-l FILE, --low_triangular_matrix FILE   Low triangular matrix of ANI values

-t FILE, --ani_table FILE               Tab separated table of ANI values

--anirb									Calculate ANI with ani.rb (should be installed separately) (slow). --input_list/--input_dir required

--fastani								Calculate ANI with fastANI (should be installed separately) (fast). --input_list/--input_dir required

--input_list FILE						List of full paths of genomes for ANI calculation

--input_dir DIRECTORY					Path (may be relative) to directory containig genomes for ANI calculation

--extension STRING						Fasta files extension for use with --input_dir, e.g. fna (default), fa, fasta

-p STRING, --prefix STRING              Prefix for output files

-m MODE, --mode MODE                    Tree inference method: UPGMA (default), NJ, both or none

-H, --heatmap                           Draw a heatmap

-A, --ascii_tree                        Draw ASCII tree to stdout

-d, --plot_dendrogram                   Plot a dendrogram

--reroot								Reroot tree at midpoint


### Input files

**Low triangular matrix** - matrix produced by [fastANI](https://github.com/ParBLiSS/FastANI) with "--matrix" option.

**ANI table** - tab-separated file of such structure: 

`Genome1    Genome2    ANI[   ...]`

May be produced by e.g. [ani.rb](https://github.com/lmrodriguezr/enveomics):

`for i in *fna; do for j in *fna; do echo -ne "$i\t$j\t"; ani.rb -q -a -1 $i -2 $j 2>/dev/null; done; done >> ani.rb.tsv`


Or by [fastANI](https://github.com/ParBLiSS/FastANI).


If you have **ani.rb** or **fastANI** installed in Your environment you may use --anirb or --fastani to calculate ANI for list of genomes using *--input_list* or folder containing genomes *--input_dir*. Genomes should have FASTA format, You may provide an extension using *--extension* option.


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
[fastANI](https://github.com/ParBLiSS/FastANI)
[ani.rb](https://github.com/lmrodriguezr/enveomics)


## Known bugs

AssertionError during midpoint rooting
`File "/somewhere/site-packages/dendropy/datamodel/treemodel.py", line 5076, in reroot_at_midpoint
	assert break_on_node is not None or target_edge is not None
AssertionError`

