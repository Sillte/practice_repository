import random
import math 
from collections import defaultdict
from collections import deque
from pprint import pprint
from itertools import zip_longest

def generate_tree(n):
    """ Generate ``Tree``. 
    """
    nodes = [None] * n
    nodes[0] = None
    nodes[1] = 0
    for i in range(1, n):
        nodes[i] = random.randint(0, i-1)
    return nodes

def generate_tree2(n):
    """ Generate ``Tree``, which is difficult to solve for Naive algorithm when n is large.  
    """

    nodes = [None] * n
    root = random.randint(0, n-1)
    nodes[root] = None
    for i in reversed(range(root)):
        nodes[i] = i+1

    for i in reversed(range(root+1, n)):
        nodes[i] = i-1
    return nodes

def generate_query(n):
    return sorted((random.randint(0, n-1), random.randint(0, n-1)))

def get_root(parents):
    return parents.index(None)

def get_ordered_nodes(parents):
    """ Return the nodes so that if i < j, then depth[i] <= depth[j] should be guaranteed.
    """
    node_to_child = defaultdict(list)
    for i, p in enumerate(parents):
        if p is not None:
            node_to_child[p].append(i)
    root = get_root(parents)
    result = [None] * len(parents)
    queue = deque([root])
    index = 0
    while queue:
        node = queue.popleft()
        result[index] = node
        index += 1
        for child in node_to_child[node]:
            queue.append(child)
    assert len(result) == index
    return result


def get_level(parents):
    """ Get the levels of tree.
    """
    levels = [None] * len(parents)
    indices = get_ordered_nodes(parents)
    levels[indices[0]] = 0
    for i in indices[1:]:
        levels[i] = levels[parents[i]] + 1
    return levels



class LCASolver:
    """ Solve Lowest Common Ancestor
    """
    def __init__(self, nodes):
        self.preprocess(nodes)

    def preprocess(self):
        """ Required Preprocess. 
        """
        raise NotImplementedError()

    def query(self, child1, child2): 
        """ Return the largest index of the node which is the common ancestor of ``child1`` and ``child2``.
        """
        raise NotImplementedError()

class Naive(LCASolver):
    """ Naive Solver
    """
    def preprocess(self, nodes):
        self.nodes = nodes
        self.levels = get_level(nodes)

    def query(self, child1, child2): 
        x, y = child1, child2
        while True: 
            if x == y:
                return x
            if self.levels[x] < self.levels[y]:
                y = self.nodes[y]
            else:
                x = self.nodes[x]

class Root(LCASolver):
    """ <O(N), O(sqrt(N))> Algorithm
    """
    def preprocess(self, parents):
        self.parents = parents
        levels = get_level(parents)

        n_level = max(levels)
        rn = math.ceil(math.sqrt(n_level))
        print("rn" ,rn)

        ps_root = [None] * len(parents)
        ordered_nodes = get_ordered_nodes(parents)
        ps_root[ordered_nodes[0]] = None
        for n in ordered_nodes[1:]:
            if levels[n] % rn != 0:
                ps_root[n] = ps_root[parents[n]]
            else:
                ps_root[n] = parents[n]
        self.ps_root = ps_root
        self.levels = levels

    def query(self, child1, child2): 
        x, y = child1, child2
        ps_root = self.ps_root
        levels = self.levels

        while ps_root[x] != ps_root[y]:
            if self.levels[x] < self.levels[y]:
                y = self.ps_root[y]
            else:
                x = self.ps_root[x]

        while True: 
            if x == y:
                return x
            if self.levels[x] < self.levels[y]:
                y = self.parents[y]
            else:
                x = self.parents[x]

class Divide(LCASolver):
    """ <O(NlogN), O(logN)> algorithm.
    """
    def preprocess(self, nodes):
        levels = get_level(nodes)

        ps_root = dict()
        """ Another implementation
        node_to_child = defaultdict(list)
        for i, p in enumerate(nodes):
            if p is not None:
                node_to_child[p].append(i)
        def _make_ps_root(node, ancesters):
            index = 0
            while (1 << index) <= len(ancesters):
                d = 1 << index 
                ps_root[node, index] = ancesters[-d]
                index += 1
            ancesters.append(node)
            for child in node_to_child[node]:
                _make_ps_root(child, ancesters)
            ancesters.pop()
        _make_ps_root(get_root(nodes), [])
        """

        max_level = max(levels)
        for i in range(len(nodes)): 
            ps_root[i, 0]  = nodes[i]

        s_max = math.floor(math.log2(max_level))
        for j in range(1, s_max + 1):
            for i in range(len(nodes)):
                if (i, j-1) in ps_root:
                    pi = ps_root[i, j-1]
                    if (pi, j-1) in ps_root:
                        ps_root[i, j]  = ps_root[pi, j-1] 
    
        self.nodes = nodes
        self.levels = levels
        self.ps_root = ps_root

    def query(self, child1, child2): 
        levels = self.levels
        ps_root = self.ps_root
        c1, c2 = sorted((child1, child2), key=lambda k: self.levels[k])
        while levels[c1] < levels[c2]:
            d = levels[c2] - levels[c1]
            c2 = ps_root[c2, math.floor(math.log2(d))]

        assert levels[c1] == levels[c2]

        while c1 != c2: 
            L = levels[c1]
            s = math.floor(math.log2(L))

            if self.nodes[c1] == self.nodes[c2]:
                c1, c2 = self.nodes[c1], self.nodes[c2]
                continue

            while 0 <= s:
                if ps_root[c1, s] != ps_root[c2, s]:
                    c1, c2 = ps_root[c1, s], ps_root[c2, s]
                    break
                s -= 1
        return c1

def check_smalltree():
    n = 1000
    parents = generate_tree(n)
    n_query = 1000
    solvers = [Naive(parents), Root(parents), Divide(parents)]
    #solvers = [Naive(parents), Root(parents)]
    for _ in range(n_query):
        node1, node2 = generate_query(len(parents))
        results = [solver.query(node1, node2) for solver in solvers] 
        print(results)
        assert all(result == results[0] for result in results), "Algorithm is not compatible."

def check_largetree():
    n = 100000
    parents = generate_tree2(n)
    print("generated", len(parents))
    n_query = 100000
    solvers = [Root(parents), Divide(parents) ]
    #solvers = [Root(parents), Naive(parents)]
    #solvers = [Convert(parents)]
    for _ in range(n_query):
        node1, node2 = generate_query(len(parents))
        results = [solver.query(node1, node2) for solver in solvers] 
        print(results)
        assert all(result == results[0] for result in results), "Algorithm is not compatible."


if __name__ == "__main__":
    #check_smalltree()
    check_largetree()
    nodes = generate_tree2(5)
    print(nodes)
    #print(nodes)
    #check_largetree()
    print(get_root(nodes))

