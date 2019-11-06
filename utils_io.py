# function to read the fna file. Returns the ids and the dictionary
# of relationships of id -> sequence for fast lookups
def read_fasta_file(filename):
    dictionary = {}
    ids = []
    with open(filename) as f:
        content = f.read().splitlines()
    currId = None
    for index, line in enumerate(content):
        if index % 2 == 0:
            currId = line[1:]
            ids.append(currId)
        else:
            dictionary[currId] = line
    return ids, dictionary
    
# writes the difference matrix to file
def write_distance_file(ids, matrix):
    count = len(ids)
    with open('genetic_distances.txt', 'w') as f:
        f.write('\t' + '\t'.join(ids) + '\n')
        for i in xrange(count):
            f.write(ids[i] + '\t' + '\t'.join(map(str, matrix[i])) + '\n')


# uses preorder traversal to traverse the tree and generate the edges file
def write_edge_file(root):
    traversal = []
    # recursive preorder traversal
    def preorder_traversal(node):
        if not node or not node.children:
            return
        for child, distance in node.children.items():
            # traverse root first
            traversal.append((node.id, child.id, distance))
            # then traverse children recursively
            preorder_traversal(child)
    preorder_traversal(root)
    # writes the traversal list to file
    with open('edges.txt', 'w') as f:
        for (parent, child, distance) in traversal:
            f.write(parent + '\t' + child + '\t' + str(distance) + '\n')


# uses post order traversal to write the newick file
def write_newick_file(original_ids, root):
    # recursive postorder traversal
    def postorder_traversal(node):
        if not node.children:
            # if it a leaf node then return the original id
            if 1 <= int(node.id) <= 61:
                return original_ids[int(node.id) - 1]
        vals = []
        # recursively traverse the children first
        for child, distance in node.children.items():
            vals.append(postorder_traversal(child) + ':' + str(distance))
        result = '(' + ','.join(vals) + ')'
        return result
    # get the root nodes children.
    # split it up into a 2 pair and 1 extra node for the root
    (first, fd), (second, sd), (third, td) = root.children.items()
    # recursively call the postorder traversal function to generate the newick format
    first_part = '(' + postorder_traversal(second) + ':' + str(sd) + ','+ postorder_traversal(third) + ':' + str(td) + ')'
    newick =  '(' + first_part + ':' + str(fd) + ');'
    with open('tree.txt', 'w') as f:
        f.write(newick)


# write the percentages in the bootstrap file
def write_bootstrap(percentages):
    with open('bootstrap.txt', 'w') as f:
        for percent in percentages:
            f.write(str(percent) + '\n')