from node import Node
import sys
from neighbour_joining import nei_saitou
import utils
from random import randint
import bootstrap

def main(filename):
    # First step is to read the file and get the sequences
    ids, sequences = utils.read_fasta_file(filename)

    # Then generate the distance matrix
    distMatrix = utils.get_distMatrix(ids, sequences)

    # write the distances.txt file (distance matrix)
    utils.write_distMatrix(ids, distMatrix)
    
    #This is used to number the sequences/nodes in the tree to be constructed
    seqCounter = 120

    root = nei_saitou(ids, distMatrix, seqCounter)

    # Writing the edge file
    utils.write_edge_file(root)

    # Writing the newick file
    utils.write_newick_file(ids, root)

    #bootstrap calculations
    #Bootstrap_instance = bootstrap.Bootstrap()
    #percentages = Bootstrap_instance.bootstrap(root, ids, sequences)
    #utils.write_bootstrap(percentages)


if __name__ == '__main__':
    main(sys.argv[1])
    