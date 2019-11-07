### Neighbor joining

Nei-Saitou [neighbor-joining algorithm](https://en.wikipedia.org/wiki/Neighbor_joining) for phylogeny construction.

### Prereq:

- Python 2.7
- R, Packages - “ape” and “RColorBrewer”
- `install.packages(c('ape','RColorBrewer'))`

### How to Run?

To run the code:
- `python main.py hw3.fna`

This generates the following files:
1. genetic_distances.txt
2. edges.txt
3. tree.txt
4. bootstrap.txt

All the above files are in the folder `submission`.

### Visualize

For visualization:

- `Rscript scripts/hw3-plot-newick.r tree.txt hw3-tip-labels.txt`
- `Rscript scripts/hw3-plot-edges.r edges.txt hw3-tip-labels.txt`

For bonus visualization:

- `Rscript scripts/hw3-plot-edges.r edges.txt hw3-tip-labels.txt boot.txt`
