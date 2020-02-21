## `ani_distance_viz.py`

Phylogenetic tree inference and heatmap drawing from ANI-derived genomic distances.

### Usage

`ani_distance_viz.py [-h] [-l LOW_TRIANGULAR_MATRIX] [-t ANI_TABLE] [-p PREFIX] [-m MODE] [-H] [-A] [-d]`

-h, --help                              Show this help message and exit

-l FILE, --low_triangular_matrix FILE   Low triangular matrix of ANI values

-t FILE, --ani_table FILE               Tab separated table of ANI values

-p STRING, --prefix STRING              Prefix for output files

-m MODE, --mode MODE                    Tree inference method: UPGMA (default), NJ, both or none

-H, --heatmap                           Draw a heatmap

-A, --ascii_tree                        Draw ASCII tree to stdout

-d, --plot_dendrogram                   Plot a dendrogram

#### Input files

Low triangular matrix - matrix produced by [fastANI](https://github.com/ParBLiSS/FastANI) with "--matrix" option.

ANI table - tab-separated file of such structure: `Genome1    Genome2    ANI[   ...]`

May be produced by e.g. [ani.rb](https://github.com/lmrodriguezr/enveomics):
`for i in *fna; do for j in *fna; do echo -ne "$i\t$j\t"; ani.rb -q -a -1 $i -2 $j 2>/dev/null; done; done >> ani.rb.tsv`


Or by [fastANI](https://github.com/ParBLiSS/FastANI).

