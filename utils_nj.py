## This file contains the util functions used by the nei-saitou neighbor-joining algorithm

# Calculates the Q matrix from the distance matrix and returns the two closest nodes and their distance
def get_qmatrix(dist_matrix):
    N = len(dist_matrix)
    qmatrix = [[0] * N for _ in range(N)]
    # keep track of the minimum value and indices in the Q matrix
    min_i = min_j = 0
    min_val = float('inf')
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            qmatrix[i][j] = (N - 2) * dist_matrix[i][j] - sum(dist_matrix[i]) - sum(dist_matrix[j])
            if qmatrix[i][j] < min_val:
                min_i = i
                min_j = j
                min_val = qmatrix[i][j]
    return min_i, min_j, min_val

# Calculates the edge lengths to the u node and returns the lengths.
def calculate_edge_lengths(dist_matrix, min_i, min_j, N):
    # Distances to the new internal node
    edge_i = 1/ 2.0 * dist_matrix[min_i][min_j] + 1 / (2.0  * (N - 2))* (sum(dist_matrix[min_i]) - sum(dist_matrix[min_j]))
    # distance from min_j node to u
    edge_j = dist_matrix[min_i][min_j] - edge_i
    return edge_i, edge_j

def get_new_dist_matrix(dist_matrix, mini, minj, N):
    updated_distances = [[0] * (N + 1) for _ in range(N + 1)]
    for i in xrange(N):
        for j in xrange(N):
            updated_distances[i][j] = dist_matrix[i][j]
    # update the distances to the new node
    for k in range(N):
        updated_distances[N][k] = (0.5) * (dist_matrix[mini][k] + dist_matrix[minj][k] - dist_matrix[mini][minj])
        updated_distances[k][N] = updated_distances[N][k]

    # Create a new distance matrix 
    new_dist_matrix = [[0] * (N - 1) for _ in range(N - 1)]
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
            new_dist_matrix[keep_i][keep_j] = updated_distances[i][j]
            keep_j += 1
        keep_i += 1

    return new_dist_matrix