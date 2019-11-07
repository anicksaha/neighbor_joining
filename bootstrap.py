from node import Node
import sys
import neighbour_joining
import utils
from random import randint

class Bootstrap:
    # do the 100 bootstrap sampling
    def bootstrap(self, root, ids, sequences):
        # original order of the ids in the original tree 
        ordered_ids = self.dfs(root,None)
        original_partitions = self.bfs(root)
        # Since we have 59 internal nodes.
        partitions_match = [0] * 59
        #Starting the 100 inferences
        for _ in range(100):
            bootstrap_sequences = {}
            indices = []
            # Using a random seq to get the length of the sequences. Assuming all the sequences to be of the same length.
            indices_to_generate = len(sequences['152801'])
            for i in range(indices_to_generate):
                indices.append(randint(0, indices_to_generate-1))

            for id in ids:
                bootstrap_seq = ''
                for index in indices:
                    bootstrap_seq += sequences[id][index]
                bootstrap_sequences[id] = bootstrap_seq

            # Neisaitou and Distance Matrix for the new sequences
            distMatrix = utils.get_distMatrix(ids, bootstrap_sequences)
            seqCounter = 120
            Neighbour_Joining_instance = neighbour_joining.Neighbour_Joining()
            root = Neighbour_Joining_instance.nei_saitou(ids, distMatrix, seqCounter)
            # get the dictionary of partitions and ids under those partitions
            partitions = self.bfs(root)

            # Comparing the original and random partitions.
            for index, id in enumerate(ordered_ids):
                if original_partitions[id] == partitions[id]:
                    partitions_match[index] += 1
        percentages = [count / 100.0 for count in partitions_match]
        return percentages


    # BFS implementation to return a partitions map when a root is given
    def bfs(self, root):
        visited = []
        queue = [root]
        partitions = {}
        while queue:
            node = queue.pop(0)
            if node.children is None:
                continue
            leaves = self.dfs(node, None)
            partitions[node.id] = leaves
            for child in node.children.keys():
                queue.append(child)

        return partitions


    # DFS to get the order of the original ids
    def dfs(self, root, visited):
        if visited is None:
            visited = []
        if root is None or root.children is None:
            return 
        visited.append(root.id)
        for child in root.children.keys():
            self.dfs(child, visited)
        return visited






