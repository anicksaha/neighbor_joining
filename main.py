import utils
from neighbor_joining import nei_saitou
import sys

def main(filename):
    # read the fna file and return the ids and id sequence mappings
    ids, id_sequences = utils.read_fasta_file(filename)

    # generate the difference matrix
    distance_matrix = utils.get_difference_matrix(ids, id_sequences)

    # write the distances.txt file (distance matrix)
    utils.write_distance_file(ids, distance_matrix)

    # run the nei saitu algorithm. Returns the root of the tree
    root = nei_saitou(ids, distance_matrix)

    # write the edge file using preorder traversal
    utils.write_edge_file(root)

    # write the newick file using postorder traversal
    utils.write_newick_file(ids, root)

    # bootstrap bonus points
    # percentages = utils.bootstrap(root, ids, id_sequences)
    # utils.write_bootstrap(percentages)

if __name__ == '__main__':
    main(sys.argv[1])