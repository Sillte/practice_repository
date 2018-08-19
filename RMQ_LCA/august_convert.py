""" Examples of conversion between RMQ and LCA.
""" 
import random
import math
from collections import deque
from collections import defaultdict
from itertools import zip_longest

def generate_sequence(n):
    """ Generate ``list``. 
    """
    return [random.randint(-n, +n) for _ in range(n)]


def generate_tree(n):
    """ Generate ``Tree``. 
    """
    nodes = [None] * n
    nodes[0] = None
    nodes[1] = 0
    for i in range(1, n):
        nodes[i] = random.randint(0, i-1)
    return nodes


def generate_query(n):
    return sorted((random.randint(0, n-1), random.randint(0, n-1)))


def get_root(parents):
    return parents.index(None)

def get_ordered_nodes(parents):
    """ Return the ``nodes`` so that if i < j, then depth[i] <= depth[j] should be guaranteed.
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

class RMQRoot:
    """ RMQSolver (Root)
    """
    def __init__(self, seq):
        self.preprocess(seq)

    def _value_min(self, args):
        return min(args, key=lambda e: (self.seq[e], e))

    def preprocess(self, seq):
        n = len(seq)
        rn = math.ceil(math.sqrt(n))
        sub_division = list(zip_longest(*([iter(range(len(seq)))] * rn), fillvalue=None))
        sub_division[-1] = [s for s in sub_division[-1] if s is not None]  # Be careful when s is 0. 
        self.seq = seq
        self.index_min = lambda i: self.seq
        self.b_size = rn
        self.b_mins = list(map(self._value_min, sub_division))

    def query(self, left, right): 
        left, right = sorted((left, right))  # For this problem setting, order adjustment is required.
        left_b = math.ceil(left / self.b_size)
        right_b = math.floor(right / self.b_size)
        if left_b + 1 < right_b:
            cand1 = range(left, self.b_size * (left_b + 1))
            cand2 = self.b_mins[left_b + 1:right_b]
            cand3 = range(right_b * self.b_size, right + 1)
            return self._value_min(map(self._value_min, (cand1, cand2, cand3)))
        else:
            return self._value_min(range(left, right + 1))


class LCARoot:
    """ <O(N), O(sqrt(N))> Algorithm
    """
    def __init__(self, parents):
        self.preprocess(parents)

    def preprocess(self, parents):
        self.parents = parents
        levels = get_level(parents)

        n_level = max(levels)
        rn = math.ceil(math.sqrt(n_level))

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

class RMQConvert:
    def __init__(self, seq):
        self.preprocess(seq)

    def preprocess(self, seq):
        tree = [None] * len(seq)
        stack = deque([])
        for index, elem in enumerate(seq): 
            if not stack:
                stack.append((index, elem))
            else:
                tmp = deque([])
                while stack and elem < stack[-1][1]: 
                    tmp.append(stack.pop())
                parent = index
                while tmp:
                    t_i, _ = tmp.pop()
                    tree[t_i] = parent
                    parent = t_i
                stack.append((index, elem))

        tmp = deque([])
        while stack:
            tmp.append(stack.pop())
        parent = None
        while tmp:
            t_i, _ = tmp.pop()
            tree[t_i] = parent
            parent = t_i
        self.tree = tree
        self.lca_solver = LCARoot(self.tree)

    def query(self, left, right):
        return self.lca_solver.query(left, right)


class LCAConvert:
    def __init__(self, parents):
        self.preprocess(parents)

    def preprocess(self, nodes):
        levels = get_level(nodes)
        node_to_child = defaultdict(list)
        for i, p in enumerate(nodes):
            if p is not None:
                node_to_child[p].append(i)

        euler_tour = list()
        """ Recursive function may failed when depth is large. 
        def _construct_tour(node):
            euler_tour.append(node)
            if node not in node_to_child:
                return
            for child in node_to_child[node]:
                _construct_tour(child)
                euler_tour.append(node)
            return
        _construct_tour(get_root(nodes))
        """
        expands_set = set()
        root = get_root(nodes)
        route = [root]
        expands_set.add(root)
        stack = deque([])
        for child in node_to_child[root]:
            stack.append(root)
            stack.append(child)
        while stack:
            current = stack.pop()
            route.append(current)
            if current in expands_set:
                continue
            expands_set.add(current)
            if not node_to_child[current]:
                continue
            for child in node_to_child[current]:
                stack.append(current)
                stack.append(child)
        self.euler_tour = route
 
        self.euler_levels = list(map(lambda t: levels[t], self.euler_tour))

        index_to_node = dict()
        for index, node in enumerate(self.euler_tour):
            if index not in index_to_node:
                index_to_node[index] = node

        self.index_to_node = index_to_node
        self.node_to_index = {value:key for key, value in index_to_node.items()}
        self.rmq_solver = RMQRoot(self.euler_levels)

    def query(self, child1, child2):
        left, right = self.node_to_index[child1], self.node_to_index[child2]
        result_index = self.rmq_solver.query(left, right) 
        return self.index_to_node[result_index]


def check_smallseq():
    n = 10000
    seq = generate_sequence(n)
    n_query = 1000
    solvers = [RMQConvert(seq), RMQRoot(seq)]
    for _ in range(n_query):
        left, right = generate_query(len(seq))
        print("left", "right", left, right)
        results = [solver.query(left,right) for solver in solvers] 
        print("results", results)
        assert all(result == results[0] for result in results), "Algorithm is not compatible."

def check_smalltree():
    n = 10000
    seq = generate_tree(n)
    n_query = 1000
    solvers = [LCAConvert(seq), LCARoot(seq)]
    for _ in range(n_query):
        left, right = generate_query(len(seq))
        print("left", "right", left, right)
        results = [solver.query(left,right) for solver in solvers] 
        print("results", results)
        assert all(result == results[0] for result in results), "Algorithm is not compatible."



if __name__ == "__main__":
    check_smallseq()
    check_smalltree()
    #seq = [2, 4, 3, 1, 6, 7, 8, 9, 1, 7]
    #solver = RMQConvert(seq)



