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