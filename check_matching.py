#!/usr/bin/python
# coding=utf8

from coxeter_graph import *
from simplicial_complex import SimplicialComplex
from matching_generator import MatchingGenerator

import os
import sys
import pickle

MATCHINGS_DIR = 'matchings'


if __name__ == '__main__':
    
    if len(sys.argv) < 3 or sys.argv[1] == "help":
        print "Usage: python %s A|B|D|E|F|H|tA|tB|tC|tD|tE|tF|tG|tI n [d] [-v|-vv] [-l]" % sys.argv[0]
        sys.exit()
    
    type = sys.argv[1]
    assert type in ['A', 'B', 'D', 'E', 'F', 'H', 'tA', 'tB', 'tC', 'tD', 'tE', 'tF', 'tG', 'tI']
    
    n = int(sys.argv[2])
    
    try:
        d = int(sys.argv[3])
    except Exception:
        d = None
    
    print "type: %s" % type
    print "n=%d" % n
    
    verbosity = 0
    if '-v' in sys.argv:
        verbosity = 1
    if '-vv' in sys.argv:
        verbosity = 2
    
    generator = MatchingGenerator(debug=False)
    graph = None
    
    try:
        if type == 'A':
            graph = SphericalACoxeterGraph(n)
        elif type == 'B':
            graph = SphericalBCoxeterGraph(n)
        elif type == 'D':
            graph = SphericalDCoxeterGraph(n)
        elif type[0] != 't':
            # spherical exceptional cases
            graph = SphericalExceptionalCoxeterGraph(type, n)
        elif type == 'tA':
            graph = AffineACoxeterGraph(n)
        elif type == 'tB':
            graph = AffineBCoxeterGraph(n)
        elif type == 'tC':
            graph = AffineCCoxeterGraph(n)
        elif type == 'tD':
            graph = AffineDCoxeterGraph(n)
        else:
            # affine exceptional cases
            graph = AffineExceptionalCoxeterGraph(type, n)
    
    except AssertionError:
        print "Unknown Coxeter type %s_%d." % (type, n)
        sys.exit()
    
    complex = SimplicialComplex(graph)
    d_values = complex.relevant_d_values() if d is None else [d]
    ranks = {}
    
    for d in d_values:
        print "*** d=%d ***" % d
        
        complex.clear_matching()
        
        try:
            generator.generate_matching(complex, d)
            matching = generator.matching
        
        except NotImplementedError:
            # MatchingGenerator does not implement this matching
            # try to load matching from file
            filename = os.path.join(MATCHINGS_DIR, "%s_%d_%d.p" % (type, n, d))
            if os.path.isfile(filename):
                # load matching from file
                with open(filename, 'r') as f:
                    matching = pickle.load(f)
            
            elif d not in complex.relevant_d_values():
                # a trivial matching works
                v = complex.vertices[0]
                matching = [(sigma, tuple(u for u in sigma if u!=v)) for sigma in complex.simplices if v in sigma]
            
            else:
                print "Matching not found."
                sys.exit()
        
        if verbosity >= 2:
            print "Matching:"
        
        for (sigma, tau) in matching:
            if verbosity >= 2:
                print sigma, tau
            complex.add_to_matching(sigma, tau, d)
        
        complex.apply_matching()
        complex.compute_morse_complex()
        complex.describe_matching(d, verbosity=verbosity)
        ranks[d] = complex.get_ranks()
    
    if '-l' in sys.argv:
        print
        print "Homology:"
        # print latex description of homology
        homology = {i: [] for i in xrange(n)}
        for d in d_values:
            for (i, rank) in enumerate(ranks[d]):
                if rank > 0:
                    homology[i].append((d, rank))
        for i in xrange(n):
            # print i-th homology group
            hom = homology[i]
            if len(hom) == 0:
                print 0
            else:
                print " \oplus ".join(("\{%d\}" % d) + ("" if rank == 1 else "^{%d}" % rank) for d, rank in hom)

