# This file contains the util functions to be used

# Function that reads the given fasta file and gets the sequences.
def read_fasta_file(filename):
    sequences = {}
    ids = []
    with open(filename) as f:
        lines = f.read().splitlines()
    curr = None
    for index, line in enumerate(lines):
        if index % 2 == 0:
            curr = line[1:]
            ids.append(curr)
        else:
            sequences[curr] = line
    return ids, sequences

# Function to calculate the distance between two sequences as a helper 
def distance(seq1,seq2):
    total = len(seq1)
    mismatch = 0
    for i in range(total):
        if seq1[i] != seq2[i]:
            mismatch += 1
    if mismatch == 0:
        return 0
    return mismatch / float(total)


# Fucntion to get the overall distance matrix for all the sequences.
def get_distMatrix(ids, sequences):
    total = len(ids)
    matrix = [[0] * total for _ in range(total)]
    for i, id1 in enumerate(ids):
        for j, id2 in enumerate(ids):
            matrix[i][j] = distance(sequences[id1], sequences[id2])
    return matrix


# writes the distance matrix to a file
def write_distMatrix(ids, matrix):
    num_ids = len(ids)
    with open('genetic_distances.txt', 'w') as f:
    	# First append the row of ids.
        f.write('\t' + '\t'.join(ids) + '\n')
        for i in range(num_ids):
            f.write(ids[i] + '\t' + '\t'.join(map(str, matrix[i])) + '\n')


# Calculates the Q matrix from the distance matrix and returns the two closest tips and their distance
def construct_QMatrix(distMatrix):
    N = len(distMatrix)
    # Matrix to store the Q values
    QMatrix = [[0] * N for _ in range(N)]

    # We also want to return the two closest tips
    mini = minj = 0
    minVal = float('inf')

    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            QMatrix[i][j] = (N - 2) * distMatrix[i][j] - sum(distMatrix[i]) - sum(distMatrix[j])
            # keep track of the minimum value and indices in the Q matrix
            if QMatrix[i][j] < minVal:
                mini = i
                minj = j
                minVal = QMatrix[i][j]
    return mini, minj, minVal


# uses post order traversal to write the newick file
def write_newick_file(seqIds, root):
    # recursive postorder traversal
    def postorder_traversal(node):
        if not node.children:
            if 1 <= int(node.id) <= 61:
                return seqIds[int(node.id) - 1]
        visited = []
        # recursively traverse the children first
        for child, distance in node.children.items():
            visited.append(postorder_traversal(child) + ':' + str(distance))
        result = '(' + ','.join(visited) + ')'
        return result
 
    (first, f), (second, s), (third, t) = root.children.items()
    first_part = '(' + postorder_traversal(second) + ':' + str(s) + ','+ postorder_traversal(third) + ':' + str(t) + ')'
    newick =  '(' + first_part + ':' + str(f) + ');'
    with open('tree.txt', 'w') as f:
        f.write(newick)


# Calculates the edge lengths to the u node and returns the lengths.
def calculate_edge_lengths(distMatrix, mini, minj, N):
    # Distances to the new internal node
    edge_i = 1/ 2.0 * distMatrix[mini][minj] + 1 / (2.0  * (N - 2))* (sum(distMatrix[mini]) - sum(distMatrix[minj]))
    # distance from minj node to u
    edge_j = distMatrix[mini][minj] - edge_i
    return edge_i, edge_j


def calculate_new_distMatrix(distMatrix, mini, minj, N):
    updatedDistances = [[0] * (N + 1) for _ in range(N + 1)]
    for i in xrange(N):
        for j in xrange(N):
            updatedDistances[i][j] = distMatrix[i][j]
    # update the distances to the new node
    for k in range(N):
        updatedDistances[N][k] = (0.5) * (distMatrix[mini][k] + distMatrix[minj][k] - distMatrix[mini][minj])
        updatedDistances[k][N] = updatedDistances[N][k]

    # Creates a new distance matrix 
    new_distMatrix = [[0] * (N - 1) for _ in range(N - 1)]
    keep_i = keep_j = 0
    for i in range(N + 1):
        # Replacing these two with a new node
        if i == mini or i == minj:
            continue
        keep_j = 0
        for j in range(N + 1):
            # Replacing these two with the new node
            if j == mini or j == minj:
                continue
            new_distMatrix[keep_i][keep_j] = updatedDistances[i][j]
            keep_j += 1
        keep_i += 1

    return new_distMatrix


# write the percentages in the bootstrap file
def write_bootstrap(percentages):
    with open('boot.txt', 'w') as f:
        for percent in percentages:
            f.write(str(percent) + '\n')

def preOrder(root, visited):
    if visited is None:
        visited = []
    if root is None or root.children is None:
        return 
    for child, distance in root.children.items():
        visited.append((root.id, child.id, distance))
        preOrder(child, visited)
    return visited


# uses preorder traversal to traverse the tree and generate the edges file
def write_edge_file(root):
    visited = preOrder(root, None)
    # Write to the edges file
    with open('edges.txt', 'w') as f:
        for (parent, child, distance) in visited:
            f.write(parent + '\t' + child + '\t' + str(distance) + '\n')

