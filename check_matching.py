#!/usr/bin/python
# coding=utf8

from coxeter_graph import *
from simplicial_complex import SimplicialComplex
from matching_generator import MatchingGenerator

import sys



if __name__ == '__main__':
    
    if len(sys.argv) < 4 or sys.argv[1] == "help":
        print "Usage: python %s [A | B | D | tA | tB | tC | tD] n d" % sys.argv[0]
        sys.exit()
    
    type = sys.argv[1]
    assert type in ['A', 'B', 'D', 'tA', 'tB', 'tC', 'tD']
    
    n = int(sys.argv[2])
    d = int(sys.argv[3])
    
    print "type: %s" % type
    print "n=%d" % n
    print "d=%d" % d
    
    generator = MatchingGenerator(debug=True)
    graph = None
    
    if type == 'A':
        graph = SphericalACoxeterGraph(n)
    
    elif type == 'B':
        graph = SphericalBCoxeterGraph(n)
    
    elif type == 'D':
        graph = SphericalDCoxeterGraph(n)
    
    elif type == 'tA':
        graph = AffineACoxeterGraph(n)
    
    elif type == 'tB':
        graph = AffineBCoxeterGraph(n)
    
    elif type == 'tC':
        graph = AffineCCoxeterGraph(n)
    
    elif type == 'tD':
        graph = AffineDCoxeterGraph(n)
    
    else:
        raise Exception("Unknown Coxeter type: %s" % type)
    
    
    complex = SimplicialComplex(graph)
    generator.generate_matching(complex, d)
    
    complex.clear_matching()

    for (sigma, tau) in generator.matching:
        # print sigma, tau
        complex.add_to_matching(sigma, tau, d)
    
    complex.apply_matching()
    complex.compute_morse_complex()
    complex.describe_matching(d)


