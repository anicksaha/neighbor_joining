import random
from neighbor_joining import nei_saitou

# Calculate the difference between two input id_sequences
# simply counts mismatches and return the percentage as a float
def calculate_difference(s1,s2):
    totalCount = len(s1)
    mismatchCount = 0
    for i in xrange(totalCount):
        if s1[i] != s2[i]:
            mismatchCount += 1
    if mismatchCount == 0:
        return 0
    return mismatchCount / float(totalCount)

# function to calculate the difference matrix from ids and the id sequence mapping
# uses the calculate_difference function
# returns the matrix
def get_difference_matrix(ids, id_sequences):
    count = len(ids)
    matrix = [[0] * count for _ in xrange(count)]
    for i, id1 in enumerate(ids):
        for j, id2 in enumerate(ids):
            matrix[i][j] = calculate_difference(id_sequences[id1], id_sequences[id2])
    return matrix


# get the leaves under this node using dfs
# returns a set of node.ids
def get_leaves(node):
    leaves = set()
    def dfs(node):
        global count
        if not node or not node.children:
            leaves.add(node.id)
        else:
            for child in node.children.keys():
                dfs(child)
    dfs(node)
    return leaves


# return the dictionary of id -> set of leaves under that node
# using bfs with a queue
def get_partitions(root):
    queue = [root]
    partition = {}
    while len(queue):
        tmp = []
        for node in queue:
            # this is a leaf
            if not node.children:
                continue
            leaves = get_leaves(node)
            partition[node.id] = leaves
            for child in node.children.keys():
                tmp.append(child)
        queue = tmp
    return partition


# get the dfs order of the tree
def get_id_order(root):
    order = []
    def dfs(node):
        if not node or not node.children:
            return
        else:
            order.append(node.id)
            for child in node.children.keys():
                dfs(child)
    dfs(root)
    return order

# do the 100 bootstrap sampling
def bootstrap(original_root, ids, id_sequences):
    # original order of teh ids in the original tree (dfs)
    original_order = get_id_order(original_root)
    original_partitions = get_partitions(original_root)
    partition_count = [0] * 59
    for _ in xrange(100):
        # new bootstrap sequence
        bootstrap_sequences = {}
        # the columns to swap
        all_indices = [i for i in xrange(len(id_sequences['152801']))]
        indices = [random.choice(all_indices) for _ in xrange(len(all_indices))]
        for id in ids:
            # choose randomly the to reconstruct the sequence
            new_sequence = ''
            for index in indices:
                new_sequence += id_sequences[id][index]
            bootstrap_sequences[id] = new_sequence
        # do the nei_saitou and get a new tree
        distance_matrix = get_difference_matrix(ids, bootstrap_sequences)
        root = nei_saitou(ids, distance_matrix)
        # get the dictionary of partitions and ids under those partitions
        partitions = get_partitions(root)

        # compare partitions to see which trees make the same partition as the original
        for index, id in enumerate(original_order):
            if original_partitions[id] == partitions[id]:
                partition_count[index] += 1
    percentages = [count / 100.0 for count in partition_count]
    return percentages