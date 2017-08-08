#!/usr/bin/python2
from random import shuffle, randint
from math import log, ceil

"""
This file contains an implementation of Munro-Patersons algorithm for exakt rank selection using multiple passes
and limited storage. 

Further optimizations could be added to the implementation:
    * Merging from multiple samplesm reduces the memory consumption by a constant factor
    * Keeping only every second element from intermediate samples reduces the memory consumption by roughly 1/2.
    * Increasing the parameter s during the iterations.
"""


def even_merge(a, b):
    """Performs the even merge operation required for Munro-Paterson.
    
    It takes two lists of sorted in ascending order and creates a new sorted list that contains every second 
    element from input lists. If both samples are i-samples, this function yields an i+1-sample.

    :param a: The first sorted list.
    :param b: The second sorted list.
    :return: A new sorted list that contains every second element from  and b. 
    """
    assert len(a) == len(b)
    assert is_sorted(a)
    assert is_sorted(b)
    m = []

    aix = 1
    bix = 1
    while aix < len(a) and bix < len(b):
        if a[aix] <= b[bix]:
            m.append(a[aix])
            aix += 2
        else:
            m.append(b[bix])
            bix += 2

    while aix < len(a):
        m.append(a[aix])
        aix += 2

    while bix < len(b):
        m.append(b[bix])
        bix += 2

    assert is_sorted(m)
    return m


def is_sorted(a):
    """ Checks whether a given list is sorted (ascending order).

    :param a: The list.
    :return: True for sorted arrays, otherwise False.
    """
    return a == sorted(a)


class ISample:
    """ This class constructs an i-sample in a streaming fashion. 

    A new object needs to be created for every new pass over a data stream.
    """

    def __init__(self, s):
        """ Creates a new object to construct an i-sample.

        :param s: The number of elements in the working sample. Parameter has the same name in the 
            Munro-Pattersen paper.
        """
        self.s = s
        # Working storage
        self.working = []

        # Subsamples
        self.buffers = []
        self.max = None
        self.count = 0

    def inspect(self, v):
        """ Processes a new element v from a stream.

        :param v: The element. 
        """
        self.working.append(v)
        if self.max is not None:
            self.max = max(self.max, v)
        else:
            self.max = v
        self.count += 1

        # If the subsample is full, we need to merge buffers
        if len(self.working) >= self.s:
            self.append_to_level(0, sorted(self.working))
            self.working = []
            self.merge()

    def append_to_level(self, level, buf):
        """ Function that registers a full i-sample in the collection. 
        
        This function should never be used by applications.

        :param level: Level of the sample to be registered.
        :param buf: The i-sample.
        """
        if len(self.buffers) <= level:
            assert (level - len(self.buffers)) == 0
            self.buffers.append([buf])
        else:
            self.buffers[level].append(buf)

    def merge(self):
        """ Merges as many i-samples as possible.
        
        This function should never be used by applications.
        """
        for level in range(0, len(self.buffers)):
            assert len(self.buffers[level]) <= 2
            if len(self.buffers[level]) == 2:
                a = self.buffers[level].pop()
                b = self.buffers[level].pop()
                self.append_to_level(level + 1, even_merge(a, b))

    def finalize(self):
        """ Fills the sample with the maximum observed value until the size reaches 2**i * s and returns the sample.
        
        :return: The i-sample, the level i
        """
        if self.count == 0:
            return [], 0

        upper_i = int(ceil(log(self.count / float(self.s), 2)))
        upper_i = max(upper_i, 0)
        target_count = 2 ** upper_i * self.s
        while self.count != target_count:
            self.inspect(self.max)

        return self.buffers[-1][0], len(self.buffers) - 1


class MunroPatersonExactMultiPass:
    """ A class that implements Munro and Patersons algorithm for exact rank computation under limited storage using
    multipe passes.
    """

    def __init__(self, target_rank, s):
        """ Creates an object to perform rank selection in limited storage.

        :param target_rank: The rank of the element, we are looking for.
        :param s: The parameter s of the algorithm, which is the length of the constructed i-samples. 
            Selecting this as m/log(n) allows us to process $n$ elements in the stream with a working storage of 
            S elements when S is sufficiently large.
            The working storage is the maximum number of elements in the lists used during the 
            construction of the i-samples.
        """
        self.target_rank = target_rank
        self.s = s

        self.upper = float("inf")
        self.u_rank = 0
        self.lower = float("-inf")
        self.l_rank = 0

        self.i_s = ISample(s)

    def inspect(self, v):
        """ Processes a value in the stream. 
        
        finalize needs to be executed after each iteration.

        :param v: The value.
        """
        if v <= self.upper:
            self.u_rank += 1
        if v <= self.lower:
            self.l_rank += 1
            return
        if v >= self.upper:
            return

        self.i_s.inspect(v)

    def finalize(self):
        """ End the current iteration and checks whether the algorithm has converged yet.
        
        :return: Converged (True / False), Bounds (exclusive lower and upper value) and the current sample level.
        """
        sample, r = self.i_s.finalize()

        # We got lucky and hit the target with the filter. Return.
        if self.u_rank == self.target_rank:
            return True, (self.upper, self.upper), r
        if self.l_rank == self.target_rank:
            return True, (self.lower, self.lower), r

        # Are we converged yet?
        if r == 0:
            val = sample[self.target_rank - self.l_rank - 1]
            return True, (val, val), r

        # We need to be very carefull with the target rank here. It is relative to the lower bound.
        rank = self.target_rank - self.l_rank

        lix = int(ceil(rank / 2.0 ** r)) - r - 1
        # If the lower bound is smaller than zero, we have a problem and just reuse the old lower filter
        if lix >= 0:
            self.lower = sample[lix]
        self.l_rank = 0

        uix = int(ceil(rank / 2.0 ** r)) - 1
        # If the upper bound is smaller than zero, we might as well use the first element from our sample
        if uix <= 0:
            uix = 0
        self.upper = sample[uix]
        self.u_rank = 0

        self.i_s = ISample(self.s)
        return False, (self.lower, self.upper), r


def test(s, i):
    """ Run a test for the Munro-Paterson algorithm with a shuffled range of integers of size 2**i * s.

    :param s: The parameter s, which also controls the size of the working sample.
    :param i: The power of two variable i.
    """
    converged = False
    tr = randint(1, 2 ** i * s + 1)
    mp = MunroPatersonExactMultiPass(tr, s)
    bounds = None
    while not converged:
        # vals = list(range(1,2**i * s + 1))
        vals = list(range(1, 2 ** i * s + 1))
        shuffle(vals)

        for v in vals:
            mp.inspect(v)
        converged, bounds, level = mp.finalize()
        print "Rank: %s, Converged: %s, Bounds: %s, Sample Level: %s" % (tr, converged, bounds, level)
    assert bounds[0] == tr


if __name__ == '__main__':
    while True:
        print "-------------"
        test(64, 10)
