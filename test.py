#!/usr/bin/python
# coding=utf8

from matching_generator import MatchingGenerator, powerset
from coxeter_graph import *
from coxeter_type import *

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
    
    


if __name__ == '__main__':
    unittest.main()


