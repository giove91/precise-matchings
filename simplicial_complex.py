#!/usr/bin/python
# coding=utf8

from itertools import combinations, chain, permutations
import copy

from complex import Cell, Edge, Complex
from coxeter_graph import CoxeterGraph

from nzmath import matrix


class Simplex:
    def __init__(self, simplicial_complex, coxeter_graph, vertices):
        self.simplicial_complex = simplicial_complex
        self.vertices = vertices
        self.weight = coxeter_graph.weight(vertices)
        
        self.up_arcs = []   # list of couples (arc, simplex)
        self.down_arcs = [] # list of couples (arc, simplex)
        
        self.matching_simplex = None
    
    
    def is_matched(self):
        return self.matching is not None
    
    def dimension(self):
        return len(self.vertices)
    
    def __str__(self):
        return u'<Simplex %s, weight %s>' % (self.vertices.__str__(), self.weight.__str__())


class SimplicialComplex:
    def __init__(self, coxeter_graph):
        self.coxeter_graph = coxeter_graph
        self.size = coxeter_graph.size
        
        if self.coxeter_graph.category == CoxeterGraph.SPHERICAL:
            # vertices are numbered 1,2,...,n
            self.vertices = range(1, self.size+1)
            # all simplices are present
            dimension = self.size
        elif self.coxeter_graph.category == CoxeterGraph.AFFINE:
            # vertices are numbered 0,1,...,n
            self.vertices = range(self.size)
            # top-dimensional simplices are not present
            dimension = self.size-1
        else:
            raise Exception("Unknown graph category")
        
        self.simplices = {sigma: Simplex(self, coxeter_graph, sigma) for sigma in chain.from_iterable(combinations(self.vertices,r) for r in xrange(dimension+1))}
        
        # Create cells
        self.cells = {vertices: Cell(simplex.dimension(), label=vertices) for (vertices, simplex) in self.simplices.iteritems()}
        
        # Create edges
        self.edges = {}
        for (vertices, simplex) in self.simplices.iteritems():
            for (i,v) in enumerate(vertices):
                # Remove vertex v
                vertices2 = tuple([w for w in vertices if w != v])
                simplex2 = self.simplices[vertices2]
                
                edge = Edge(self.cells[vertices], self.cells[vertices2], (-1)**i)
                self.edges[vertices, vertices2] = edge
        
        # Create a complex
        self.complex = Complex(self.cells.values(), self.edges.values())
        
        self.matching = set()
        self.morse_complex = None
        self.is_matching_applied = False
    
    
    def clear_matching(self):
        for s in self.simplices.itervalues():
            s.matching_simplex = None
        
        self.matching = set()
        self.is_matching_applied = False
        for e in self.complex.edges:
            if e.is_in_matching:
                e.remove_from_matching()
    
    
    def add_to_matching(self, sigma, tau, d):
        sigma = tuple(sorted(sigma))
        tau = tuple(sorted(tau))
        
        if self.simplices[sigma].weight.component(d) != self.simplices[tau].weight.component(d):
            raise Exception("Trying to match simplices with different weight: %s weight %d and %s weight %d" % (str(self.simplices[sigma].vertices), self.simplices[sigma].weight.component(d), str(self.simplices[tau].vertices), self.simplices[tau].weight.component(d)))
        
        for x in permutations([self.simplices[sigma], self.simplices[tau]]):
            assert x[0].matching_simplex is None
            x[0].matching_simplex = x[1]
        
        self.matching.add((sigma,tau))
    
    
    def import_matching_from_complex(self, d):
        for e in self.complex.edges:
            if e.is_in_matching:
                self.add_to_matching(e.high.label, e.low.label, d)
        self.is_matching_applied = True
    
    
    def apply_matching(self, debug=False):
        """
        Apply matching to self.complex.
        If the matching is not acyclic, an exception is raised.
        """
        for e in self.complex.edges:
            if (e.high.label, e.low.label) in self.matching:
                # add cell to matching
                e.add_to_matching()
                if not self.complex.is_acyclic(e.low, e.high.d, print_cycle=debug):
                    if debug:
                        print e.high, e.low
                    raise Exception("Matching is not acyclic")
        self.is_matching_applied = True
    
    
    def compute_morse_complex(self):
        """
        Compute the Morse complex.
        """
        if not self.is_matching_applied:
            self.apply_matching()
        self.morse_complex = self.complex.morse_reduction()
    
    
    def is_matching_precise(self, d, debug=False):
        if self.morse_complex is None:
            self.compute_morse_complex()
        
        for e in self.morse_complex.edges:
            assert e.deg != 0 # only edges with incidence != 0 are stored
            
            s1 = self.simplices[e.high.label]
            s2 = self.simplices[e.low.label]
            
            # check if these simplices are relevant
            if all(self.coxeter_graph.is_simplex_relevant(s.vertices) for s in [s1, s2]):
                if debug:
                    print s1, s2, d
                # check if their weights differ by 1
                if s1.weight.component(d) != s2.weight.component(d) + 1:
                    return False
        return True
    
    
    def critical_simplices(self):
        for s in self.simplices.itervalues():
            if s.matching_simplex is None:
                yield s
    
    
    def get_ranks(self):
        """
        Returns the ranks over Q of the boundaries of the Morse complex, from the 1-dim to the top-dim.
        """
        assert self.morse_complex is not None
        r, boundaries = self.morse_complex.get_boundaries()
        ranks = []

        for k in sorted(boundaries.iterkeys()):
            boundary = boundaries[k]
            b = boundary.copy() if (boundary.row > 0 and boundary.column > 0) else matrix.Matrix(row=1, column=1)
            b.toFieldMatrix()
            ranks.append(b.rank())
        
        return ranks
        
    
    def relevant_d_values(self):
        """
        Return list of relevant values for d.
        """
        relevant_d = set()
        for simplex in self.simplices.itervalues():
            for (d, w) in simplex.weight.w.iteritems():
                if w >= 1:
                    relevant_d.add(d)
        return list(sorted(relevant_d))
    
    
    def describe_matching(self, d, verbosity=0):
        if verbosity >= 1:
            print "Critical simplices:"
            for s in self.critical_simplices():
                print s.vertices, "\t", "w=%d" % s.weight.component(d)
        
        if self.morse_complex is None:
            self.compute_morse_complex()
        
        if verbosity >= 2:
            print "In the Morse complex there are %d edge(s):" % len(self.morse_complex.edges)
            for e in self.morse_complex.edges:
                print e
        
        if self.is_matching_precise(d):
            print "The matching is precise."
            ranks = self.get_ranks()
            print "Ranks (from 1-dim to %d-dim):" % len(ranks), ranks
        else:
            print "The matching is *not* precise."

        


