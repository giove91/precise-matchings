#!/usr/bin/python
# coding=utf8

from matching_generator import MatchingGenerator, powerset
from coxeter_graph import *
from coxeter_type import *
from complex import Cell, Edge, Complex
from simplicial_complex import SimplicialComplex

import unittest
import fractions

def phi(n):
    return sum(1 for k in xrange(1, n+1) if fractions.gcd(n, k) == 1)


class TestCoxeterType(unittest.TestCase):
    
    def test_weight(self):
        w = Weight({2: 3, 3: 2, 4: 1, 5: 1, 6: 1})
        self.assertEqual(w.component(2), 3)
        self.assertEqual(w.component(3), 2)
        self.assertEqual(w.component(4), 1)
        self.assertEqual(w.component(5), 1)
        self.assertEqual(w.component(6), 1)
        self.assertEqual(w.component(7), 0)
        self.assertEqual(w.component(8), 0)
    
    
    def test_coxeter_type(self):
        t = CoxeterType('A', 5)
        self.assertEqual(t.weight(), Weight({2: 3, 3: 2, 4: 1, 5: 1, 6: 1}))

        u = CoxeterType('A', 3)
        self.assertEqual(u.weight(), Weight({2: 2, 3: 1, 4: 1}))
        
        self.assertEqual(t.weight() + u.weight(), Weight({2: 5, 3: 3, 4: 2, 5: 1, 6: 1}))
    
    
    def test_coxeter_type_A(self):
        for n in xrange(1, 40):
            t = CoxeterType('A', n)
            degree = n*(n+1)/2
            self.assertEqual(sum(phi(d)*w for (d,w) in t.weight().w.iteritems()), degree)
    
    def test_coxeter_type_B(self):
        for n in xrange(2, 40):
            t = CoxeterType('B', n)
            degree = n*(n+1) - n
            self.assertEqual(sum(phi(d)*w for (d,w) in t.weight().w.iteritems()), degree)
    
    def test_coxeter_type_D(self):
        for n in xrange(4, 40):
            t = CoxeterType('D', n)
            degree = n*(n-1)
            self.assertEqual(sum(phi(d)*w for (d,w) in t.weight().w.iteritems()), degree)
    
    def test_coxeter_type_E8(self):
        t = CoxeterType('E', 8)
        degree = 1+7+11+13+17+19+23+29  # sum of exponents
        self.assertEqual(sum(phi(d)*w for (d,w) in t.weight().w.iteritems()), degree)
    
    def test_exceptional_coxeter_type(self):
        m = 44
        exponents = {
                ('E', 6): [1,4,5,7,8,11],
                ('E', 7): [1,5,7,9,11,13,17],
                ('E', 8): [1,7,11,13,17,19,23,29],
                ('F', 4): [1,5,7,11],
                ('H', 3): [1,5,9],
                ('H', 4): [1,11,19,29],
                ('I', 2): [1, m-1]
            }
        
        for (type, n) in exponents.iterkeys():
            t = CoxeterType(type, n, m if type == 'I' else None)
            degree = sum(exponents[type, n])
            self.assertEqual(sum(phi(d)*w for (d,w) in t.weight().w.iteritems()), degree)



class TestCoxeterGraph(unittest.TestCase):
    
    def test_create_graph_from_arcs(self):
        graph = CoxeterGraph()
        graph.size = 5
        arcs = [(2,3,3), (3,4,4), (2,4,5), (8,9,3)]
        graph.create_graph_from_arcs(arcs)
        
        self.assertEqual(len(graph.vertices), 5)
        self.assertEqual(len(graph.arcs), 4)
    
    
    def test_find_connected_components(self):
        graph = CoxeterGraph()
        graph.size = 5
        arcs = [(2,3,3), (3,4,4), (2,4,5), (8,9,3)]
        graph.create_graph_from_arcs(arcs)
        
        self.assertEqual(len(graph.vertices), 5)
        self.assertEqual(len(graph.arcs), 4)
        
        components = graph.find_connected_components((2,3,4,8,9))
        self.assertEqual(components, [(8,9), (2,3,4)])
        
        components = graph.find_connected_components((2,3,4,9))
        self.assertEqual(components, [(9,), (2,3,4)])
        
        components = graph.find_connected_components((2,4,9))
        self.assertEqual(components, [(9,), (2,4)])
        
        components = graph.find_connected_components((2,3,4))
        self.assertEqual(components, [(2,3,4)])
        
        components = graph.find_connected_components((2,4))
        self.assertEqual(components, [(2,4)])
        
        components = graph.find_connected_components((2,))
        self.assertEqual(components, [(2,)])
        
        components = graph.find_connected_components(())
        self.assertEqual(components, [])
    
    
    def test_spherical_A_graph(self):
        graph = SphericalACoxeterGraph(5)
        self.assertEqual(len(graph.vertices), 5)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [1,2,3,4,5])
        self.assertEqual(len(graph.arcs), 4)
        self.assertEqual(graph.category, CoxeterGraph.SPHERICAL)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (1,3,4,5)
        self.assertEqual(graph.weight(s), Weight({2: 3, 3: 1, 4: 1}))

        s = (1,2,4,5)
        self.assertEqual(graph.weight(s), Weight({2: 2, 3: 2}))

        s = (2,3,4)
        self.assertEqual(graph.weight(s), Weight({2: 2, 3: 1, 4: 1}))
    
    
    def test_spherical_B_graph(self):
        graph = SphericalBCoxeterGraph(5)
        self.assertEqual(len(graph.vertices), 5)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [1,2,3,4,5])
        self.assertEqual(len(graph.arcs), 4)
        self.assertEqual(graph.category, CoxeterGraph.SPHERICAL)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (1,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 3).weight())

        s = (1,2,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('B', 2).weight() + CoxeterType('A', 2).weight())

        s = (2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight())
        
        s = (1,2,3)
        self.assertEqual(graph.weight(s), CoxeterType('B', 3).weight())
    
    
    def test_spherical_D_graph(self):
        graph = SphericalDCoxeterGraph(6)
        self.assertEqual(len(graph.vertices), 6)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [1,2,3,4,5,6])
        self.assertEqual(len(graph.arcs), 5)
        self.assertEqual(graph.category, CoxeterGraph.SPHERICAL)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (1,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 4).weight())

        s = (1,2,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 1).weight() + CoxeterType('A', 2).weight())

        s = (2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight())
        
        s = (1,2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('D', 4).weight())
        
        s = (1,2,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('D', 5).weight())
        
        s = (1,2,3,4,5,6)
        self.assertEqual(graph.weight(s), CoxeterType('D', 6).weight())
    
    
    def test_spherical_E8_graph(self):
        graph = SphericalExceptionalCoxeterGraph('E', 8)
        self.assertEqual(len(graph.vertices), 8)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [1,2,3,4,5,6,7,8])
        self.assertEqual(len(graph.arcs), 7)
        self.assertEqual(graph.category, CoxeterGraph.SPHERICAL)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (1,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 3).weight())

        s = (1,2,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 2).weight() + CoxeterType('A', 2).weight())

        s = (2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight())
        
        s = (1,2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 4).weight())
        
        s = (1,2,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 5).weight())
        
        s = (1,2,3,4,5,6)
        self.assertEqual(graph.weight(s), CoxeterType('A', 6).weight())
        
        s = (1,2,5,6,8)
        self.assertEqual(graph.weight(s), CoxeterType('A', 2).weight() + CoxeterType('A', 3).weight())
        
        s = (1,4,5,6,7,8)
        self.assertEqual(graph.weight(s), CoxeterType('D', 5).weight() + CoxeterType('A', 1).weight())
        
        s = (1,3,4,5,6,7,8)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('E', 6).weight())
        
        s = (2,3,4,5,6,7,8)
        self.assertEqual(graph.weight(s), CoxeterType('E', 7).weight())
        
        s = (1,2,3,4,5,6,7,8)
        self.assertEqual(graph.weight(s), CoxeterType('E', 8).weight())
        
    
    def test_spherical_F4_graph(self):
        graph = SphericalExceptionalCoxeterGraph('F', 4)
        self.assertEqual(len(graph.vertices), 4)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [1,2,3,4])
        self.assertEqual(len(graph.arcs), 3)
        self.assertEqual(graph.category, CoxeterGraph.SPHERICAL)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (1,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 2).weight())

        s = (1,2,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 2).weight())

        s = (2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('B', 3).weight())
        
        s = (1,2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('F', 4).weight())
        
        s = (1,2,3)
        self.assertEqual(graph.weight(s), CoxeterType('B', 3).weight())
    
    
    def test_spherical_H4_graph(self):
        graph = SphericalExceptionalCoxeterGraph('H', 4)
        self.assertEqual(len(graph.vertices), 4)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [1,2,3,4])
        self.assertEqual(len(graph.arcs), 3)
        self.assertEqual(graph.category, CoxeterGraph.SPHERICAL)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (1,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 2).weight())

        s = (1,2,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('I', 2, 5).weight())

        s = (2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight())
        
        s = (1,2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('H', 4).weight())
        
        s = (1,2,3)
        self.assertEqual(graph.weight(s), CoxeterType('H', 3).weight())
    
    
    def test_spherical_I2_graph(self):
        m = 102
        graph = SphericalExceptionalCoxeterGraph('I', 2, m)
        self.assertEqual(len(graph.vertices), 2)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [1,2])
        self.assertEqual(len(graph.arcs), 1)
        self.assertEqual(graph.category, CoxeterGraph.SPHERICAL)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (1,)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight())

        s = (2,)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight())

        s = (1, 2)
        self.assertEqual(graph.weight(s), CoxeterType('I', 2, m).weight())
    
    
    def test_affine_A_graph(self):
        graph = AffineACoxeterGraph(5)
        self.assertEqual(len(graph.vertices), 6)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [0,1,2,3,4,5])
        self.assertEqual(len(graph.arcs), 6)
        self.assertEqual(graph.category, CoxeterGraph.AFFINE)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (1,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 3).weight())

        s = (1,2,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 2).weight() + CoxeterType('A', 2).weight())

        s = (2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight())
        
        s = (0,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 2).weight())
        
        s = (0,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight())
    
    
    def test_affine_B_graph(self):
        graph = AffineBCoxeterGraph(5)
        self.assertEqual(len(graph.vertices), 6)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [0,1,2,3,4,5])
        self.assertEqual(len(graph.arcs), 5)
        self.assertEqual(graph.category, CoxeterGraph.AFFINE)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (1,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 3).weight())

        s = (1,2,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 2).weight() + CoxeterType('A', 1).weight() + CoxeterType('A', 1).weight())

        s = (2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight())
        
        s = (0,1)
        self.assertEqual(graph.weight(s), CoxeterType('B', 2).weight())
        
        s = (0,1,2)
        self.assertEqual(graph.weight(s), CoxeterType('B', 3).weight())
        
        s = (1,2,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('D', 5).weight())
        
        s = (0,)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight())
    
    
    def test_affine_C_graph(self):
        graph = AffineCCoxeterGraph(5)
        self.assertEqual(len(graph.vertices), 6)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [0,1,2,3,4,5])
        self.assertEqual(len(graph.arcs), 5)
        self.assertEqual(graph.category, CoxeterGraph.AFFINE)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (1,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('B', 3).weight())

        s = (1,2,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 2).weight() + CoxeterType('B', 2).weight())

        s = (2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight())
        
        s = (0,1)
        self.assertEqual(graph.weight(s), CoxeterType('B', 2).weight())
        
        s = (1,2,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('B', 5).weight())
        
        s = (0,)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight())
    
    
    def test_affine_D_graph(self):
        graph = AffineDCoxeterGraph(5)
        self.assertEqual(len(graph.vertices), 6)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [0,1,2,3,4,5])
        self.assertEqual(len(graph.arcs), 5)
        self.assertEqual(graph.category, CoxeterGraph.AFFINE)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (1,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 3).weight())

        s = (1,2,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 1).weight() + CoxeterType('A', 2).weight())

        s = (2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight())
        
        s = (0,1)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 1).weight())
        
        s = (1,2,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('D', 5).weight())
        
        s = (0,)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight())
        
        s = (0,1,2,3)
        self.assertEqual(graph.weight(s), CoxeterType('D', 4).weight())
    
    
    def test_affine_D_graph2(self):
        graph = AffineDCoxeterGraph(4)
        self.assertEqual(len(graph.vertices), 5)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [0,1,2,3,4])
        self.assertEqual(len(graph.arcs), 4)
        self.assertEqual(graph.category, CoxeterGraph.AFFINE)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (1,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 1).weight() + CoxeterType('A', 1).weight())

        s = (1,2,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight())

        s = (2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight())
        
        s = (0,1)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 1).weight())
        
        s = (1,2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('D', 4).weight())
        
        s = (0,)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight())
        
        s = (0,1,2,3)
        self.assertEqual(graph.weight(s), CoxeterType('D', 4).weight())
        
        s = (0,2,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight())
    
    
    def test_affine_E6_graph(self):
        graph = AffineExceptionalCoxeterGraph('tE', 6)
        self.assertEqual(len(graph.vertices), 7)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [0,1,2,3,4,5,6])
        self.assertEqual(len(graph.arcs), 6)
        self.assertEqual(graph.category, CoxeterGraph.AFFINE)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (1,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 1).weight() + CoxeterType('A', 2).weight())

        s = (1,2,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight() + CoxeterType('A', 1).weight())

        s = (2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight())
        
        s = (1,2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 4).weight())
        
        s = (1,2,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('D', 5).weight())
        
        s = (1,2,3,4,5,6)
        self.assertEqual(graph.weight(s), CoxeterType('E', 6).weight())
        
        s = (0,1,2,5,6)
        self.assertEqual(graph.weight(s), CoxeterType('A', 5).weight())
        
        s = (0,1,3,4,5,6)
        self.assertEqual(graph.weight(s), CoxeterType('A', 2).weight() + CoxeterType('A', 2).weight() + CoxeterType('A', 2).weight())
        
        s = (0,1,2,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('E', 6).weight())
        
        s = (1,2,3,5)
        self.assertEqual(graph.weight(s), CoxeterType('D', 4).weight())


    def test_affine_E7_graph(self):
        graph = AffineExceptionalCoxeterGraph('tE', 7)
        self.assertEqual(len(graph.vertices), 8)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [0,1,2,3,4,5,6,7])
        self.assertEqual(len(graph.arcs), 7)
        self.assertEqual(graph.category, CoxeterGraph.AFFINE)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (1,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 3).weight())

        s = (1,2,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 2).weight() + CoxeterType('A', 2).weight())

        s = (2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight())
        
        s = (1,2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 4).weight())
        
        s = (1,2,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 5).weight())
        
        s = (1,2,3,4,5,6)
        self.assertEqual(graph.weight(s), CoxeterType('A', 6).weight())
        
        s = (1,2,5,6,7)
        self.assertEqual(graph.weight(s), CoxeterType('A', 2).weight() + CoxeterType('A', 2).weight() + CoxeterType('A', 1).weight())
        
        s = (1,3,4,5,6,7)
        self.assertEqual(graph.weight(s), CoxeterType('A', 5).weight() + CoxeterType('A', 1).weight())
        
        s = (2,3,4,5,6,7)
        self.assertEqual(graph.weight(s), CoxeterType('D', 6).weight())
        
        s = (1,2,3,4,5,7)
        self.assertEqual(graph.weight(s), CoxeterType('E', 6).weight())
        
        s = (1,2,3,4,5,6,7)
        self.assertEqual(graph.weight(s), CoxeterType('E', 7).weight())
        
        s = (0,1,2,3,4,7)
        self.assertEqual(graph.weight(s), CoxeterType('D', 6).weight())


    def test_affine_E8_graph(self):
        graph = AffineExceptionalCoxeterGraph('tE', 8)
        self.assertEqual(len(graph.vertices), 9)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [0,1,2,3,4,5,6,7,8])
        self.assertEqual(len(graph.arcs), 8)
        self.assertEqual(graph.category, CoxeterGraph.AFFINE)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (1,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 3).weight())

        s = (1,2,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 2).weight() + CoxeterType('A', 2).weight())

        s = (2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight())
        
        s = (1,2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 4).weight())
        
        s = (1,2,3,4,5)
        self.assertEqual(graph.weight(s), CoxeterType('A', 5).weight())
        
        s = (1,2,3,4,5,6)
        self.assertEqual(graph.weight(s), CoxeterType('A', 6).weight())
        
        s = (1,2,5,6,8)
        self.assertEqual(graph.weight(s), CoxeterType('A', 2).weight() + CoxeterType('A', 3).weight())
        
        s = (1,4,5,6,7,8)
        self.assertEqual(graph.weight(s), CoxeterType('D', 5).weight() + CoxeterType('A', 1).weight())
        
        s = (1,3,4,5,6,7,8)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('E', 6).weight())
        
        s = (2,3,4,5,6,7,8)
        self.assertEqual(graph.weight(s), CoxeterType('E', 7).weight())
        
        s = (1,2,3,4,5,6,7,8)
        self.assertEqual(graph.weight(s), CoxeterType('E', 8).weight())
        
        s = (0,1,2,3,4,5,6,8)
        self.assertEqual(graph.weight(s), CoxeterType('D', 8).weight())
    
    
    def test_affine_F4_graph(self):
        graph = AffineExceptionalCoxeterGraph('tF', 4)
        self.assertEqual(len(graph.vertices), 5)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [0,1,2,3,4])
        self.assertEqual(len(graph.arcs), 4)
        self.assertEqual(graph.category, CoxeterGraph.AFFINE)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (1,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 2).weight())

        s = (1,2,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('B', 2).weight())

        s = (2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('A', 3).weight())
        
        s = (1,2,3,4)
        self.assertEqual(graph.weight(s), CoxeterType('B', 4).weight())
        
        s = (1,2,3)
        self.assertEqual(graph.weight(s), CoxeterType('B', 3).weight())
        
        s = (0,1,2,3)
        self.assertEqual(graph.weight(s), CoxeterType('F', 4).weight())
    
    
    def test_affine_G2_graph(self):
        graph = AffineExceptionalCoxeterGraph('tG', 2)
        self.assertEqual(len(graph.vertices), 3)
        self.assertEqual(list(sorted(graph.vertices.iterkeys())), [0,1,2])
        self.assertEqual(len(graph.arcs), 2)
        self.assertEqual(graph.category, CoxeterGraph.AFFINE)
        
        s = ()
        self.assertEqual(graph.weight(s), Weight())
        
        s = (0,1)
        self.assertEqual(graph.weight(s), CoxeterType('I', 2, 6).weight())

        s = (1,2)
        self.assertEqual(graph.weight(s), CoxeterType('A', 2).weight())

        s = (0,2)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight() + CoxeterType('A', 1).weight())
        
        s = (1,)
        self.assertEqual(graph.weight(s), CoxeterType('A', 1).weight())
    

class TestComplex(unittest.TestCase):
    
    def test_torus(self):
        # S^1 x S^1
        cells = [Cell(0), Cell(1), Cell(1), Cell(2)]
        edges = [Edge(cells[1],cells[0],0), Edge(cells[2],cells[0],0), Edge(cells[3],cells[1],0), Edge(cells[3],cells[2],0)]
        C = Complex(cells, edges)
        with self.assertRaises(AssertionError):
            edges[0].add_to_matching()
    
    def test_torus2(self):
        # S^1 x S^1
        cells = [Cell(0), Cell(1), Cell(1), Cell(1), Cell(2), Cell(2)]
        edges = [Edge(cells[4], cells[1], 1), Edge(cells[4], cells[2], 1), Edge(cells[4], cells[3], 1)] + [Edge(cells[5], cells[1], 1), Edge(cells[5], cells[2], 1), Edge(cells[5], cells[3], 1)]
        C = Complex(cells, edges)
        
        edges[5].add_to_matching()
        M = C.morse_reduction()
        self.assertEqual(M.d, 2)
        self.assertEqual(len(M.cells[0]), 1)
        self.assertEqual(len(M.cells[1]), 2)
        self.assertEqual(len(M.cells[2]), 1)
        
        ranks, boundaries = M.get_boundaries()
    
    def test_simplex(self):
        # {1,2,3}
        cells = {s: Cell(len(s)-1, label=s) for s in [(), (1,), (2,), (3,), (1,2), (1,3), (2,3), (1,2,3)]}
        edges = {(s,t): Edge(cells[s], cells[t], x) for (s,t,x) in [((1,), (), 1), ((2,), (), 1), ((3,), (), 1), ((1,2), (1,), -1), ((1,2), (2,), 1), ((1,3), (1,), -1), ((1,3), (3,), 1), ((2,3), (2,), -1), ((2,3), (3,), 1), ((1,2,3), (1,2), 1), ((1,2,3), (1,3), -1), ((1,2,3), (2,3), 1)]}
        C = Complex(cells.values(), edges.values())
        
        edges[(1,), ()].add_to_matching()
        edges[(1,2), (2,)].add_to_matching()
        edges[(1,3), (3,)].add_to_matching()
        edges[(1,2,3), (2,3)].add_to_matching()
        
        M = C.morse_reduction()
        self.assertEqual(M.d, -1)
        self.assertTrue(all(len(y) == 0 for y in M.cells.itervalues()))


class TestSimplicialComplex(unittest.TestCase):
    
    def test_with_spherical_A_graph(self):
        graph = SphericalACoxeterGraph(3)
        c = SimplicialComplex(graph)
        
        self.assertEqual(c.vertices, [1,2,3])
        self.assertEqual(set(c.simplices.iterkeys()), set([(), (1,), (2,), (3,), (1,2), (1,3), (2,3), (1,2,3)]))

        c.apply_matching()
    
    
    def test_with_affine_A_graph(self):
        graph = AffineACoxeterGraph(3)
        c = SimplicialComplex(graph)
        
        self.assertEqual(c.vertices, [0,1,2,3])
        self.assertEqual(set(c.simplices.iterkeys()), set([(), (0,), (1,), (2,), (3,), (0,1), (0,2), (0,3), (1,2), (1,3), (2,3), (0,1,2), (0,1,3), (0,2,3), (1,2,3)]))

        c.apply_matching()
    
    
    def test_matching(self):
        graph = SphericalACoxeterGraph(3)
        c = SimplicialComplex(graph)
        
        self.assertEqual(c.vertices, [1,2,3])
        self.assertEqual(set(c.simplices.iterkeys()), set([(), (1,), (2,), (3,), (1,2), (1,3), (2,3), (1,2,3)]))
        d = 3
        
        c.add_to_matching((1,), (), d)
        c.apply_matching()
        c.clear_matching()
        
        c.add_to_matching((1,), (), d)
        c.add_to_matching((1,3), (3,), d)
        c.apply_matching()
        c.clear_matching()
        
        c.add_to_matching((1,), (), d)
        c.add_to_matching((1,3), (3,), d)
        c.add_to_matching((1,2,3), (2,3), d)
        c.apply_matching()
        
        self.assertEqual(set([s.vertices for s in c.critical_simplices()]), set([(2,), (1,2)]))
        c.clear_matching()
        
        c.add_to_matching((1,), (), d)
        c.add_to_matching((1,3), (3,), d)
        c.add_to_matching((1,2,3), (2,3), d)
        
        with self.assertRaises(AssertionError):
            # a simplex is already matched
            c.add_to_matching((2,), (), d)
    
    
    def test_matching2(self):
        graph = SphericalACoxeterGraph(3)
        c = SimplicialComplex(graph)
        
        self.assertEqual(c.vertices, [1,2,3])
        self.assertEqual(set(c.simplices.iterkeys()), set([(), (1,), (2,), (3,), (1,2), (1,3), (2,3), (1,2,3)]))
        
        c.apply_matching()
        c.clear_matching()
        
        d = 5 # all simplices have weight 0
        
        c.add_to_matching((1,), (), d)
        c.apply_matching()
        c.clear_matching()
        
        c.add_to_matching((1,), (), d)
        c.add_to_matching((1,3), (3,), d)
        c.apply_matching()
        c.clear_matching()
        
        c.add_to_matching((1,), (), d)
        c.add_to_matching((1,3), (3,), d)
        c.add_to_matching((1,2,3), (2,3), d)
        c.apply_matching()
        c.clear_matching()
        
        c.add_to_matching((1,), (), d)
        c.add_to_matching((1,3), (3,), d)
        c.add_to_matching((1,2,3), (2,3), d)
        c.add_to_matching((1,2), (2,), d)
        c.apply_matching()

    
    def test_matching3(self):
        graph = SphericalACoxeterGraph(5)
        c = SimplicialComplex(graph)
        
        self.assertEqual(c.vertices, [1,2,3,4,5])
        
        c.apply_matching()
        c.clear_matching()
        
        d = 5 # all simplices have weight 0
        
        c.add_to_matching((1,), (), d)
        c.apply_matching()
        c.clear_matching()
        
        c.add_to_matching((1,), (), d)
        c.add_to_matching((1,3), (3,), d)
        c.apply_matching()
        c.clear_matching()
        
        c.add_to_matching((1,), (), d)
        c.add_to_matching((1,3), (3,), d)
        c.add_to_matching((1,2,3), (1,2), d)
        c.apply_matching()
        c.clear_matching()
        
        c.add_to_matching((1,), (), d)
        c.add_to_matching((1,3), (3,), d)
        c.add_to_matching((1,2,3), (1,2), d)
        c.add_to_matching((1,2,4), (2,4), d)
        c.apply_matching()
        c.clear_matching()
        
        c.add_to_matching((1,), (), d)
        c.add_to_matching((1,3), (3,), d)
        c.add_to_matching((1,2,3), (1,2), d)
        c.add_to_matching((1,2,4), (2,4), d)
        c.add_to_matching((2,3,4), (2,3), d) # this creates a cycle
        
        with self.assertRaises(Exception):
            c.apply_matching()
    
    
    def test_ranks(self):
        graph = SphericalACoxeterGraph(3)
        c = SimplicialComplex(graph)
        
        self.assertEqual(c.vertices, [1,2,3])
        self.assertEqual(set(c.simplices.iterkeys()), set([(), (1,), (2,), (3,), (1,2), (1,3), (2,3), (1,2,3)]))
        d = 2
        
        c.add_to_matching((1,2), (1,), d)
        c.add_to_matching((1,2,3), (1,3), d)
        c.add_to_matching((2,3), (2,), d)
        c.apply_matching()
        
        self.assertEqual(set([s.vertices for s in c.critical_simplices()]), set([(), (3,)]))
        self.assertTrue(c.is_matching_precise(d))
        
        self.assertEqual(c.get_ranks(), [1])
    
    def test_ranks2(self):
        graph = SphericalACoxeterGraph(5)
        c = SimplicialComplex(graph)
        
        self.assertEqual(c.vertices, [1,2,3,4,5])
        d = 3
        
        c.add_to_matching((1,2,3), (1,2), d)
        c.add_to_matching((1,2,3,4), (1,2,4), d)
        c.add_to_matching((1,2,3,5), (1,2,5), d)
        c.add_to_matching((1,2,3,4,5), (1,2,4,5), d)
        c.add_to_matching((1,), (), d)
        c.add_to_matching((1,3), (3,), d)
        c.add_to_matching((1,4), (4,), d)
        c.add_to_matching((1,5), (5,), d)
        c.add_to_matching((1,3,4), (3,4), d)
        c.add_to_matching((1,3,5), (3,5), d)
        c.add_to_matching((1,4,5), (4,5), d)
        c.add_to_matching((1,3,4,5), (3,4,5), d)
        c.add_to_matching((2,3,4,5), (2,3,5), d)
        c.add_to_matching((2,3,4), (2,3), d)
        c.add_to_matching((2,4), (2,), d)
        c.apply_matching()
        
        self.assertEqual(set([s.vertices for s in c.critical_simplices()]), set([(2,5), (2,4,5)]))
        self.assertTrue(c.is_matching_precise(d))
        
        self.assertEqual(c.get_ranks(), [0,0,1])


class TestMatchingGenerator(unittest.TestCase):
    
    def test_powerset(self):
        self.assertEqual(list(powerset(range(3))), [(), (0,), (1,), (2,), (0,1), (0,2), (1,2), (0,1,2)])
    
    
    def test_add_to_matching(self):
        generator = MatchingGenerator()
        generator.add_to_matching((0,2,3), 2)
        self.assertEqual(generator.matching, [((0,2,3), (0,3))])
        
        generator.add_to_matching((0,2,3), 1)
        self.assertEqual(generator.matching, [((0,2,3), (0,3)), ((0,1,2,3), (0,2,3))])
    
    
    def test_generate_spherical_A_matching(self):
        generator = MatchingGenerator()
        for n in xrange(1,9):
            for d in xrange(2,n+2):
                graph = SphericalACoxeterGraph(n)
                complex = SimplicialComplex(graph)
                generator.generate_matching(complex, d)
                
                for (sigma, tau) in generator.matching:
                    complex.add_to_matching(sigma, tau, d)
                
                self.assertEqual(sum(1 for s in complex.critical_simplices()), 2 if n%d in [0, d-1] else 0)
                self.assertTrue(complex.is_matching_precise(d))
    
    
    def test_generate_spherical_A_matching_first_present(self):
        generator = MatchingGenerator(debug=False)
        for n in xrange(0, 9):
            for d in xrange(2, 2*n+3):
                for f in xrange(n+1):
                    # print n,d,f
                    graph = SphericalACoxeterGraph(n)
                    complex = SimplicialComplex(graph)
                    generator.generate_matching(complex, d, f=f)
                    
                    for (sigma, tau) in generator.matching:
                        complex.add_to_matching(sigma, tau, d)
                    
                    simplices_not_considered = 2**n - 2**(n-f) # number of simplices not satisfying 1,2,...,f in sigma
                    # print sum(1 for s in complex.critical_simplices()) - simplices_not_considered
                    self.assertEqual(sum(1 for s in complex.critical_simplices()) - simplices_not_considered, 1 if f == n else 0 if f%d==d-1 else 2 if n%d in [f%d, d-1] else 0)
    
    
    def test_generate_spherical_A_matching_fg(self):
        generator = MatchingGenerator(debug=False)
        for n in xrange(1, 7):
            for d in xrange(2, n+2):
                for f in xrange(n+1):
                    for g in xrange(n+1):
                        # print n,d,f,g
                        graph = SphericalACoxeterGraph(n, f=f, g=g)
                        complex = SimplicialComplex(graph)
                        generator.generate_matching(complex, d, f=f, g=g)
                        
                        for (sigma, tau) in generator.matching:
                            complex.add_to_matching(sigma, tau, d)
                        
                        self.assertTrue(complex.is_matching_precise(d, debug=False))
                        
                        simplices_not_considered = 2**n - 2**(max(n-f-g,0)) # number of simplices not satisfying 1,2,...,f,n-g+1,...,n in sigma
                        num_critical = sum(1 for s in complex.critical_simplices()) - simplices_not_considered
                        
                        if g == 0:
                            self.assertEqual(num_critical, 1 if f == n else 0 if f%d==d-1 else 2 if n%d in [f%d, d-1] else 0)
    
    
    def test_generate_spherical_B_matching_odd(self):
        # case d odd
        generator = MatchingGenerator()
        for n in xrange(2,10):
            for d in xrange(3,n+2,2):
                graph = SphericalBCoxeterGraph(n)
                complex = SimplicialComplex(graph)
                generator.generate_matching(complex, d)
                
                for (sigma, tau) in generator.matching:
                    complex.add_to_matching(sigma, tau, d)
                
                self.assertEqual(sum(1 for s in complex.critical_simplices()), 0)
                self.assertTrue(complex.is_matching_precise(d))
    
    
    def test_generate_spherical_B_matching2(self):
        # case d=2
        generator = MatchingGenerator()
        d = 2
        for n in xrange(2,10):
            graph = SphericalBCoxeterGraph(n)
            complex = SimplicialComplex(graph)
            generator.generate_matching(complex, d)
            
            for (sigma, tau) in generator.matching:
                complex.add_to_matching(sigma, tau, d)
            
            self.assertEqual(sum(1 for s in complex.critical_simplices()), 2*n)
            
            self.assertTrue(complex.is_matching_precise(d))
    
    
    def test_generate_spherical_B_matching_even(self):
        # case d >= 4 even
        generator = MatchingGenerator(debug=False)
        for n in xrange(2, 10):
            for d in xrange(4, 2*n+2, 2):
                graph = SphericalBCoxeterGraph(n)
                complex = SimplicialComplex(graph)
                generator.generate_matching(complex, d)
                
                for (sigma, tau) in generator.matching:
                    complex.add_to_matching(sigma, tau, d)
                
                self.assertTrue(complex.is_matching_precise(d))
    
    
    def test_generate_spherical_B_matching_even_g(self):
        # case d >= 4 even, with g >= 1
        generator = MatchingGenerator(debug=False)
        for n in xrange(2, 9):
            for d in xrange(2, 2*n+2, 2):
                for g in xrange(1, d/2+1):
                    # print n, d, g
                    graph = SphericalBCoxeterGraph(n, g=g)
                    complex = SimplicialComplex(graph)
                    generator.generate_matching(complex, d, g=g)
                    
                    for (sigma, tau) in generator.matching:
                        complex.add_to_matching(sigma, tau, d)
                    
                    self.assertTrue(complex.is_matching_precise(d, debug=False))
    
    
    def test_generate_spherical_D_matching_odd(self):
        # case d odd
        generator = MatchingGenerator()
        for n in xrange(4,10):
            for d in xrange(3,n+2,2):
                graph = SphericalDCoxeterGraph(n)
                complex = SimplicialComplex(graph)
                generator.generate_matching(complex, d)
                
                for (sigma, tau) in generator.matching:
                    complex.add_to_matching(sigma, tau, d)
                
                self.assertEqual(sum(1 for s in complex.critical_simplices()), (2 if n%d in [0,1] else 0))
                self.assertTrue(complex.is_matching_precise(d))
    
    
    def test_generate_spherical_D_matching_even(self):
        # case d >= 6 even
        generator = MatchingGenerator(debug=False)
        for n in xrange(4,10):
            for d in xrange(6,2*n+2,2):
                # print n, d
                graph = SphericalDCoxeterGraph(n)
                complex = SimplicialComplex(graph)
                generator.generate_matching(complex, d)
                
                for (sigma, tau) in generator.matching:
                    complex.add_to_matching(sigma, tau, d)
                
                num_critical = None
                
                if n >= 5:
                    num_critical = sum((1 if k in [n,n-1] else 2) for k in xrange(4,n+1) if (n-k)%d in [0,1]) if n%d in [0, d/2+1] else 0
                    num_critical += 2 if n%d in [0,4] else 0
                    num_critical += 2 if n%d in [0,1] else 0
                else:
                    num_critical = 2 if d == 6 else 0
                
                self.assertEqual(sum(1 for s in complex.critical_simplices()), num_critical)
                self.assertTrue(complex.is_matching_precise(d))
                
    
    def test_generate_spherical_D_matching4(self):
        # case d=4
        generator = MatchingGenerator(debug=False)
        d = 4
        for n in xrange(4,10):
            # print n,d
            graph = SphericalDCoxeterGraph(n)
            complex = SimplicialComplex(graph)
            generator.generate_matching(complex, d)
            
            for (sigma, tau) in generator.matching:
                complex.add_to_matching(sigma, tau, d)
            
            num_critical = sum((1 if k in [n,n-1] else 2) for k in xrange(4,n+1) if (n-k)%d in [0,1]) if n%d in [0, d/2+1] else 0
            num_critical += 2 if n%4 in [0,1] else 0
            num_critical += 2 * (1 if n == 4 else 2 if n%4 in [0,3] else 0)
            num_critical += 1 if n == 4 else 0
            
            self.assertEqual(sum(1 for s in complex.critical_simplices()), num_critical)
            self.assertTrue(complex.is_matching_precise(d, debug=False))
            
    
    
    def test_generate_spherical_D_matching2(self):
        # case d=2
        generator = MatchingGenerator(debug=False)
        d = 2
        for n in xrange(4,10):
            # print n,d
            graph = SphericalDCoxeterGraph(n)
            complex = SimplicialComplex(graph)
            generator.generate_matching(complex, d)
            
            for (sigma, tau) in generator.matching:
                complex.add_to_matching(sigma, tau, d)
            
            self.assertEqual(sum(1 for s in complex.critical_simplices()), 2 if n%2==1 else 4)
            self.assertTrue(complex.is_matching_precise(d))
    
    
    def test_generate_spherical_D_matching_even_g(self):
        # case g >= 1, d even
        generator = MatchingGenerator(debug=False)
        for n in xrange(4,8):
            for d in xrange(2,2*n+6,2):
                for g in xrange(1, min(d/2+1, n+1)):
                    # print n, d, g
                    graph = SphericalDCoxeterGraph(n, g=g)
                    complex = SimplicialComplex(graph)
                    generator.generate_matching(complex, d, g=g)
                    
                    for (sigma, tau) in generator.matching:
                        complex.add_to_matching(sigma, tau, d)
                    
                    self.assertTrue(complex.is_matching_precise(d))
    
    
    def test_generate_affine_A_matching(self):
        generator = MatchingGenerator(debug=False)
        for n in xrange(1,9):
            for d in xrange(2,n+2):
                graph = AffineACoxeterGraph(n)
                complex = SimplicialComplex(graph)
                generator.generate_matching(complex, d)
                
                for (sigma, tau) in generator.matching:
                    complex.add_to_matching(sigma, tau, d)
                
                self.assertTrue(complex.is_matching_precise(d))
                
                # check that the result coincide with the one already computed in literature
                complex.compute_morse_complex()
                ranks = complex.get_ranks()
                
                expected_ranks = [0] * len(ranks)
                for j in xrange((n+1)/d):
                    i = n - 2*j - 1 if (n+1)%d == 0 else n - 2*j - 2
                    expected_ranks[i] = d-1 if (n+1)%d == 0 else 1
                
                # print ranks, expected_ranks
                self.assertEqual(ranks, expected_ranks)

    
    def test_generate_affine_B_matching_odd(self):
        # case d odd
        generator = MatchingGenerator()
        for n in xrange(3,9):
            for d in xrange(3,n+2,2):
                graph = AffineBCoxeterGraph(n)
                complex = SimplicialComplex(graph)
                generator.generate_matching(complex, d)
                
                for (sigma, tau) in generator.matching:
                    complex.add_to_matching(sigma, tau, d)
                
                self.assertEqual(sum(1 for s in complex.critical_simplices()), 1)
                self.assertTrue(complex.is_matching_precise(d))
    
    
    def test_generate_affine_B_matching_even(self):
        # case d even
        generator = MatchingGenerator(debug=False)
        for n in xrange(3,9):
            for d in xrange(4,2*n+2,2):
                # print n,d
                graph = AffineBCoxeterGraph(n)
                complex = SimplicialComplex(graph)
                generator.generate_matching(complex, d)
                
                for (sigma, tau) in generator.matching:
                    complex.add_to_matching(sigma, tau, d)
                
                if d >= 6:
                    num_critical = 0
                    if n % d/2 in [0, 1]:
                        num_critical += 1
                    if n % d/2 == 0:
                        num_critical += 1
                    num_critical += n / (d/2)
                    
                    if n % d/2 == 0:
                        num_critical += (n+1) / (d/2)
                    else:
                        pass
                
                for s in complex.critical_simplices():
                    # check that the quantity "dimension - weight" is constant and equal to n - n/(d/2)
                    self.assertEqual(s.dimension() - graph.weight(s.vertices).component(d), n - n/(d/2))
                
                self.assertTrue(complex.is_matching_precise(d, debug=False))
    
    
    def test_generate_affine_B_matching2(self):
        # case d=2
        generator = MatchingGenerator(debug=False)
        d = 2
        for n in xrange(3,9):
            # print n,d
            graph = AffineBCoxeterGraph(n)
            complex = SimplicialComplex(graph)
            generator.generate_matching(complex, d)
            
            for (sigma, tau) in generator.matching:
                complex.add_to_matching(sigma, tau, d)
            
            self.assertEqual(sum(1 for s in complex.critical_simplices()), 3*n if n%2==1 else 3*n+1)
            self.assertTrue(complex.is_matching_precise(d))
            
            for s in complex.critical_simplices():
                # check that the quantity "dimension - weight" is constant and equal to n - n/(d/2)
                self.assertEqual(s.dimension() - graph.weight(s.vertices).component(d), n - n/(d/2))
    
    
    def test_generate_affine_C_matching_odd(self):
        # case d odd
        generator = MatchingGenerator()
        for n in xrange(3,10):
            for d in xrange(3,n+2,2):
                graph = AffineCCoxeterGraph(n)
                complex = SimplicialComplex(graph)
                generator.generate_matching(complex, d)
                
                for (sigma, tau) in generator.matching:
                    complex.add_to_matching(sigma, tau, d)
                
                self.assertEqual(sum(1 for s in complex.critical_simplices()), 1)
                self.assertTrue(complex.is_matching_precise(d))
    
    
    def test_generate_affine_C_matching_even(self):
        # case d even
        generator = MatchingGenerator()
        for n in xrange(3,10):
            for d in xrange(2,2*n+2,2):
                graph = AffineCCoxeterGraph(n)
                complex = SimplicialComplex(graph)
                generator.generate_matching(complex, d)
                
                for (sigma, tau) in generator.matching:
                    complex.add_to_matching(sigma, tau, d)
                
                self.assertTrue(complex.is_matching_precise(d))
    
    
    def test_generate_affine_C_matching2(self):
        # case d=2
        generator = MatchingGenerator(debug=False)
        d = 2
        for n in xrange(3,10):
            graph = AffineCCoxeterGraph(n)
            complex = SimplicialComplex(graph)
            generator.generate_matching(complex, d)
            
            for (sigma, tau) in generator.matching:
                complex.add_to_matching(sigma, tau, d)
            
            # print sum(1 for s in complex.critical_simplices())
            self.assertEqual(sum(1 for s in complex.critical_simplices()), n**2+n+1)
            self.assertTrue(complex.is_matching_precise(d))
    
    
    def test_generate_affine_D_matching_odd(self):
        # case d odd
        generator = MatchingGenerator()
        for n in xrange(4,10):
            for d in xrange(3,n+2,2):
                graph = AffineDCoxeterGraph(n)
                complex = SimplicialComplex(graph)
                generator.generate_matching(complex, d)
                
                for (sigma, tau) in generator.matching:
                    complex.add_to_matching(sigma, tau, d)
                
                self.assertEqual(sum(1 for s in complex.critical_simplices()), 3 if n%d in [0,1] else 1)
                self.assertTrue(complex.is_matching_precise(d))
        
    
    def test_generate_affine_D_matching2(self):
        # case d=2, n>=5
        generator = MatchingGenerator(debug=False)
        d = 2
        for n in xrange(5,10):
            # print n,d
            graph = AffineDCoxeterGraph(n)
            complex = SimplicialComplex(graph)
            generator.generate_matching(complex, d)
            
            for (sigma, tau) in generator.matching:
                complex.add_to_matching(sigma, tau, d)
            
            self.assertEqual(sum(1 for s in complex.critical_simplices()), n+7 if n%2==0 else n+2)
            self.assertTrue(complex.is_matching_precise(d))
    
    
    def test_generate_affine_D_matching4(self):
        # case d=4, n>=5
        generator = MatchingGenerator(debug=False)
        d = 4
        for n in xrange(5,10):
            # print n,d
            graph = AffineDCoxeterGraph(n)
            complex = SimplicialComplex(graph)
            generator.generate_matching(complex, d)
            
            for (sigma, tau) in generator.matching:
                complex.add_to_matching(sigma, tau, d)
            
            self.assertTrue(complex.is_matching_precise(d))
    
    
    def test_generate_affine_D_matching_even(self):
        # case d>=6 even, n>=5
        generator = MatchingGenerator(debug=False)
        for n in xrange(5,10):
            for d in xrange(6,2*n+2,2):
                # print n,d
                graph = AffineDCoxeterGraph(n)
                complex = SimplicialComplex(graph)
                generator.generate_matching(complex, d)
                
                for (sigma, tau) in generator.matching:
                    complex.add_to_matching(sigma, tau, d)
                
                self.assertTrue(complex.is_matching_precise(d, debug=False))
    
    
    def test_generate_affine_D_matching_even_all(self):
        # case d even
        generator = MatchingGenerator(debug=False)
        for n in xrange(4,9):
            for d in xrange(2,2*n+8,2):
                # print n, d
                graph = AffineDCoxeterGraph(n)
                complex = SimplicialComplex(graph)
                generator.generate_matching(complex, d)
                
                for (sigma, tau) in generator.matching:
                    complex.add_to_matching(sigma, tau, d)
                
                self.assertTrue(complex.is_matching_precise(d))



if __name__ == '__main__':
    unittest.main()


