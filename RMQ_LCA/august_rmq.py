import random
import math 
from pprint import pprint
from itertools import zip_longest

def generate_sequence(n):
    """ Generate ``list``. 
    """
    return [random.randint(-n, +n) for _ in range(n)]

def generate_query(n):
    return sorted((random.randint(0, n-1), random.randint(0, n-1)))

class RMQSolver:
    """ Solve Range Minumization Query.
    """
    def __init__(self, seq):
        self.preprocess(seq)

    def preprocess(self):
        """ Required Preprocess. 
        """
        raise NotImplementedError()

    def query(self, left, right): 
        """ Return the minimum value between ``left`` and ``right``, inclusive.
        """
        raise NotImplementedError()

class Naive(RMQSolver):
    """ <O(N^2), O(1)> algorithm.
    """
    def preprocess(self, seq):
        self.seq = seq
        n = len(seq)
        result = [[None] * n for _ in range(n)]
        for i in range(n):
            result[i][i] = i
        for i in range(n):
            for j in range(i+1, n):
                result[i][j] = min((result[i][j-1], j), key=lambda i:(self.seq[i], i))
        self.memo = result

    def query(self, left, right): 
        return self.memo[left][right]


class Root(RMQSolver):
    """ <O(N^2), O(sqrt(N))> algorithm.
    """
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
        left_b = math.ceil(left / self.b_size)
        right_b = math.floor(right / self.b_size)
        if left_b + 1 < right_b:
            cand1 = range(left, self.b_size * (left_b + 1))
            cand2 = self.b_mins[left_b + 1:right_b]
            cand3 = range(right_b * self.b_size, right + 1)
            return self._value_min(map(self._value_min, (cand1, cand2, cand3)))
        else:
            return self._value_min(range(left, right + 1))


class SparseTree(RMQSolver):
    """ <O(NlogN), O(1)> algorithm.
    """
    def _value_min(self, args):
        return min(args, key=lambda e: (self.seq[e], e))

    def preprocess(self, seq):
        self.seq = seq
        n = len(seq)
        d = math.floor(math.log2(n)) + 1
        memo = [[None] * d for _ in range(n)]

        for i in range(n):
            memo[i][0] = i
        for j in range(1, d):
            for i in range(n):
                memo[i][j] = self._value_min((memo[i][j-1], memo[min(i+2**(j-1), n-1)][j-1]))
        self.memo = memo
    
    def query(self, left, right): 
        length = right - left
        if length == 0:
            return self.memo[right][0]

        interval_index = math.floor(math.log2(length)) 
        d = 2 ** interval_index
        assert left <= left + d
        assert right - d + 1 <= left + (d - 1) + 1

        return self._value_min((self.memo[left][interval_index], self.memo[right- d + 1][interval_index]))

class SegmentTree(RMQSolver):
    """ <O(logN), O(logN)> algorithm.
    """
    def _value_min(self, args):
        return min(args, key=lambda e: (self.seq[e], e))

    def preprocess(self, seq):
        self.seq = seq

        n = len(seq)
        heap = dict() # Index starts from 1.
        h_to_range = dict()
        def _construct(h, r):
            h_to_range[h] = r
            if r[1] == r[0]:
                heap[h] = r[0]
                return r[0]
            center = (r[1] + r[0]) // 2 
            s = _construct(2 * h, (r[0], center))
            t = _construct(2 * h + 1,(center + 1, r[1]))
            heap[h] = self._value_min((s, t))
            return heap[h]
        _construct(1, (0, n - 1))
        self.heap = heap
        self.h_to_range = h_to_range

    def query(self, left, right):
        def _search(h, r):
            p, q = r
            tp, tq = self.h_to_range[h]

            if (p, q) == (tp, tq):
                return self.heap[h]

            center = (tp + tq) // 2
            assert tp <= p <= q <= tq
            if center + 1 <= p:
                return _search(2 * h + 1, r)  
            elif q <= center:
                return _search(2 * h, r)  
            else:
                cand1 = _search(2 * h + 1, (center + 1, q))
                cand2 = _search(2 * h, (p,  center))
                return self._value_min((cand1, cand2))

        return _search(1, (left, right))
                

def check_smallseq():
    n = 1000
    seq = generate_sequence(n)
    n_query = 1000
    solvers = [Naive(seq), Root(seq), SparseTree(seq), SegmentTree(seq)]
    for _ in range(n_query):
        left, right = generate_query(len(seq))
        print("left", "right", left, right)
        results = [solver.query(left,right) for solver in solvers] 
        print("results", results)
        assert all(result == results[0] for result in results), "Algorithm is not compatible."

def check_largeseq():
    n = 10000
    seq = generate_sequence(n)
    print(seq)
    n_query = 10000
    solvers = [Root(seq), SparseTree(seq), SegmentTree(seq)]
    for _ in range(n_query):
        left, right = generate_query(len(seq))
        print("left", "right", left, right)
        results = [solver.query(left,right) for solver in solvers] 
        print("results", results)
        assert all(result == results[0] for result in results), "Algorithm is not compatible."



if __name__ == "__main__":
    #check_smallseq()
    check_largeseq()
    seq = [1 ,2, 3, 4, 5]
    solver = SegmentTree(seq)
    pprint(solver.heap)
    pprint(solver.h_to_range)
    pprint(solver.query(1, 3))
