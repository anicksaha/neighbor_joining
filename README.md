To run the code:
python main.py scripts-data/hw3.fna

This generates the following files by default:
1. genetic_distances.txt
2. edges.txt
3. tree.txt
4. bootstrap.txt
All the above files are in the folder answers.

Running the bootstrap code take a little longer. To avoid that,everything that
follow the comment #bootstrap calculations# in main.py can be commented out.

For visualization:

Rscript scripts-data/hw3-plot-newick.r tree.txt scripts-data/hw3-tip-labels.txt
or
Rscript scripts-data/hw3-plot-edges.r edges.txt scripts-data/hw3-tip-labels.txt

For bonus visualization:

Rscript scripts-data/hw3-plot-edges.r edges.txt scripts-data/hw3-tip-labels.txt boot.txt
