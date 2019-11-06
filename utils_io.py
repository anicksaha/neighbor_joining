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