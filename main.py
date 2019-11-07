from node import Node
import neighbour_joining
import utils
import bootstrap
import sys

def main(filename):
    # First step is to read the file and get the sequences
    ids, sequences = utils.read_fasta_file(filename)

    # Then generate the distance matrix
    distance_matrix = utils.get_distance_matrix(ids, sequences)

    # write the distances.txt file (distance matrix)
    utils.write_distance_matrix(ids, distance_matrix)
    
    # This is used to number nodes in the tree to be constructed
    sequence_counter = 120
    # Run the nei saitu algorithm. Returns the root of the tree
    root = neighbour_joining.nei_saitou(ids, distance_matrix, sequence_counter)

    # Writing the edge file
    utils.write_edge_file(root)

    # Writing the newick file
    utils.write_newick_file(ids, root)

    # bootstrap calculations
    percentages = bootstrap.bootstrap(root, ids, sequences)
    bootstrap.write_bootstrap(percentages)

if __name__ == '__main__':
    main(sys.argv[1])