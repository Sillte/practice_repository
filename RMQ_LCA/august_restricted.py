""" RMQ for Restricted Version.
Adjacent elements differs in +1 or -1.   
This condition corresponds to the instance which is converted from LCA.  
"""
import random
import math 
from collections import defaultdict
from collections import deque
from pprint import pprint
from itertools import zip_longest

def generate_restricted_problem(n):
    """ Generate ``list``. 
    """
    start = random.randint(-n, +n)
    result = [start]
    for _ in range(n-1):
        result.append(result[-1] + random.choice([+1, -1]))
    return result

def generate_query(n):
    return sorted((random.randint(0, n-1), random.randint(0, n-1)))

class SparseTree:
    """ <O(NlogN), O(1)> algorithm.
    """
    def __init__(self, seq):
        self.preprocess(seq)

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

class RestrictedSolver:
    """ <O(N), O(1)> algorithm.
    """
    def _value_min(self, args):
        return min(args, key=lambda e: (self.seq[e], e))

    def __init__(self, seq):
        self.seq = seq
        self._miniseq_dict = dict()
        self.preprocess(seq)

    def preprocess(self, seq):
        n = len(seq)
        bin_size = math.floor(math.log2(n) / 2) + 1
        bins = list(zip_longest(*([iter(range(n))]*bin_size), fillvalue=None))
        bins[-1] = [s for s in bins[-1] if s is not None] # s may be 0.
        bin_seq = list(map(self._value_min, bins))
        # For bin_seq, Sparse Tree is applied.
        bn = len(bin_seq)
        bd = math.floor(math.log2(bn)) + 1
        bin_memo = [[None] * bd for _ in range(bn)]
        for i in range(bn):
            bin_memo[i][0] = bin_seq[i]
        for j in range(1, bd):
            for i in range(bn):
                bin_memo[i][j] = self._value_min((bin_memo[i][j-1], bin_memo[min(i+2**(j-1), bn-1)][j-1]))

        self.bin_size = bin_size 
        self.bin_memo = bin_memo
        self.bin_seq = bin_seq

        self.miniseq_map = dict()
        for index, one_bin in enumerate(bins):
            self.miniseq_map[index] = self._calc_miniseq_map(one_bin)


    def _calc_miniseq_map(self, one_bin):
        """ Memorization for the block sequence
        Notice that the answer of index is within the block. 
        So, you have to add the start of the block to this answer.

        It is important that you can store the all patterns becase 
        all adjacent numbers differ in +1 or -1.
        """ 
        if len(one_bin) == 1:
            return [[0]]

        nums = [self.seq[i] for i in one_bin]
        min_num = lambda args: min(args, key=lambda a: nums[a])
        key = tuple([nums[i] - nums[i-1] for i in range(1, len(nums))])
        if key in self._miniseq_dict:
            #print("Memorized")
            return self._miniseq_dict[key]

        result = [[None] * len(nums) for _ in range(len(nums))]
        for i in range(len(nums)):
            result[i][i] = i

        for i in range(0, len(one_bin)): 
            for j in range(i+1, len(one_bin)): 
                result[i][j] = min_num((result[i][j-1], j))

        self._miniseq_dict[key] = result
        return result

    def _bin_query(self, left, right): 
        left_bin = left // self.bin_size
        right_bin = right // self.bin_size
        assert left_bin == right_bin, "Should be same bin"
        bin_index = left_bin
        memo_dict = self.miniseq_map[bin_index]
        # Local problem.
        left, right = left % self.bin_size, right % self.bin_size
        d_result = memo_dict[left][right]
        bin_result = left_bin * self.bin_size + d_result 
        return bin_result


    def query(self, left, right):
        left_bin = left // self.bin_size
        right_bin = right // self.bin_size
        result = right

        length = right_bin - left_bin
        if right_bin - left_bin >= 2:
            if right_bin - left_bin == 2:
                cand = self.bin_seq[left_bin + 1]
                result = self._value_min((result, cand))
            else:
                r_inner = right_bin - 1
                l_inner = left_bin + 1
                interval = math.floor(math.log2(r_inner  - l_inner))
                cands = (self.bin_memo[l_inner][interval], self.bin_memo[r_inner - (2**interval) + 1][interval])
                cand = self._value_min(cands) 
                result =  self._value_min((result, cand))

            # Processing of the edge region.
            left_cand = self._bin_query(left, (left // self.bin_size + 1) * self.bin_size - 1)
            right_cand = self._bin_query((right // self.bin_size) * self.bin_size, right)
            result = self._value_min((result, left_cand, right_cand))

        elif length == 1:
            left_cand = self._bin_query(left, (left // self.bin_size + 1) * self.bin_size - 1)
            right_cand = self._bin_query((right // self.bin_size) * self.bin_size, right)
            result = self._value_min((result, left_cand, right_cand))

        elif length == 0:
            result = self._bin_query(left, right)

        return result


def check_smallseq():
    for _ in range(10):
        n = 100000
        seq = generate_restricted_problem(n)
        n_query = 1000
        solvers = [SparseTree(seq), RestrictedSolver(seq)]
        #solvers = [SparseTree(seq)]
        for _ in range(n_query):
            left, right = generate_query(len(seq))
            #print("left", "right", left, right)
            results = [solver.query(left,right) for solver in solvers] 
            #print("results", results)
            #print("problem", seq)
            assert all(result == results[0] for result in results), "Algorithm is not compatible."


if __name__ ==  "__main__":
    check_smallseq()
    problem = generate_restricted_problem(10)
    print("problem", problem)
    solver = RestrictedSolver(problem)
    print(solver.query(3, 4))
