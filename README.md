## `genomic_distance_viz.py`

Phylogenetic tree inference and heatmap drawing from ANI (average nucleotide identity) or AAI (average aminoacid identity)-derived genomic distances.


## Usage

`genomic_distance_viz.py [-h] (-l LOW_TRIANGULAR_MATRIX | -t ANI_TABLE | --anirb | --aairb | --mummer | --fastani) [--input-list INPUT_LIST | --input-dir INPUT_DIR] [--extension EXTENSION] [-p PREFIX] [-m MODE] [--threads THREADS] [-H] [-A] [-d] [--reroot]`


### Command-line options

```
-h, --help								Show this help message and exit

-l FILE, --low-triangular-matrix FILE   Low triangular matrix of ANI values

-t FILE, --table FILE					Tab separated table of similarity (ANI or AAI) between genomes

--anirb									Calculate ANI with ani.rb (should be installed separately) (slow). --input-list/--input-dir required

--fastani								Calculate ANI with fastANI (should be installed separately) (fast). --input_list/--input_dir required

--mummer								Calculate ANI with mummer. --input_list/--input_dir required

--input-list FILE						List of full paths of genomes for ANI calculation

--input-dir DIRECTORY					Path (may be relative) to directory containig genomes for ANI calculation

--extension STRING						Fasta files extension for use with --input_dir, e.g. fna (default), fa, fasta

-p STRING, --prefix STRING				Prefix for output files

-m MODE, --mode MODE					Tree inference method: UPGMA (default), NJ, both or none

--threads THREADS						Number of CPU threads (where possible)

-H, --heatmap							Draw a heatmap

-A, --ascii-tree						Draw ASCII tree to stdout

-d, --plot-dendrogram					Plot a dendrogram

--reroot								Reroot tree at midpoint. May cause errors or incorrect trees
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
[fastANI](https://github.com/ParBLiSS/FastANI)
[ani.rb](https://github.com/lmrodriguezr/enveomics)
[aai.rb](https://github.com/lmrodriguezr/enveomics)
[Mummer4](https://github.com/mummer4/mummer)

## Known bugs

AssertionError during midpoint rooting
`File "/somewhere/site-packages/dendropy/datamodel/treemodel.py", line 5076, in reroot_at_midpoint
	assert break_on_node is not None or target_edge is not None
AssertionError`

