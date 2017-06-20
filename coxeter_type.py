#!/usr/bin/python
# coding=utf8

import copy


class Weight:
    def __init__(self, w={}):
        self.w = {d: e for (d,e) in w.iteritems() if e>0} # dictionary: d-component => exponent
    
    def component(self, d):
        """
        Returns the exponent of the d-th cyclotomic polynomial.
        """
        if d in self.w:
            return self.w[d]
        else:
            return 0
    
    def __str__(self):
        return self.w.__str__()
    
    def __repr__(self):
        return self.w.__repr__()
    
    def __eq__(self, other):
        return self.w == other.w
    
    def __neq__(self, other):
        return not self == other
    
    def __add__(self, other):
        s = set(self.w.keys()) | set(other.w.keys())
        return Weight({d: self.component(d) + other.component(d) for d in s})
    
    def __radd__(self, other):
        if other == 0:
            return copy.deepcopy(self)
        else:
            return self + other



class CoxeterType:
    """
    Type of a finite Coxeter graph.
    """
    
    def __init__(self, type, n=None, m=None):
        self.type = type
        assert self.type in ['A', 'B', 'D', 'E', 'F', 'H', 'I']
        self.n = n
        self.m = m # Relevant only in the I_2(m) case
    
    def __eq__(self, other):
        return (self.type, self.n, self.m) == (other.type, other.n, other.m)
    
    def weight(self):
        """
        Returns the weight (Poincaré polynomial of the Coxeter group,
        expressed as a product of cyclotomic polynomials).
        """
        n = self.n
        
        if self.type == 'A':
            return Weight({d: (n+1)/d for d in xrange(2, n+2)})
        
        elif self.type == 'B':
            return Weight({d: n/d if d%2 == 1 else 2*n/d for d in xrange(2, 2*n+1)})
        
        elif self.type == 'D':
            return Weight({d: n/d if d%2 == 1 else (2*(n-1)/d + (1 if n%d==0 else 0)) for d in xrange(2, 2*n-1)})
        
        else:
            exponents = {
                ('E', 6): [1,4,5,7,8,11],
                ('E', 7): [1,5,7,9,11,13,17],
                ('E', 8): [1,7,11,13,17,19,23,29],
                ('F', 4): [1,5,7,11],
                ('H', 3): [1,5,9],
                ('H', 4): [1,11,19,29],
                ('I', 2): [1, self.m-1 if self.m is not None else None]
            }[self.type, self.n]
            
            return sum(Weight({d: 1 if (e+1)%d == 0 else 0 for d in xrange(2,e+2)}) for e in exponents)


