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
def get_distance(seq1,seq2):
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
            matrix[i][j] = get_distance(sequences[id1], sequences[id2])
    return matrix


# writes the distance matrix to a file
def write_distMatrix(ids, matrix):
    num_ids = len(ids)
    with open('genetic_distances.txt', 'w') as f:
    	# First append the row of ids.
        f.write('\t' + '\t'.join(ids) + '\n')
        for i in range(num_ids):
            f.write(ids[i] + '\t' + '\t'.join(map(str, matrix[i])) + '\n')

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