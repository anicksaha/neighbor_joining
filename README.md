### Neighbor joining

Nei-Saitou [neighbor-joining algorithm](https://en.wikipedia.org/wiki/Neighbor_joining) for phylogeny construction.

### Prereq:

- Python 2.7
- R, packages in R: {“ape”, “RColorBrewer”}

`install.packages(c('ape','RColorBrewer'))`

### How to Run?

- `python main.py <path to sequence file>`
- `python main.py hw3.fna`

### Visualize

- `Rscript hw3-plot-newick.r tree.txt hw3-tip-labels.txt`
- `Rscript hw3-plot-edges.r edges.txt hw3-tip-labels.txt`
