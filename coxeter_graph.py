#!/usr/bin/python
# coding=utf8

from coxeter_type import CoxeterType, Weight


class Vertex:
    def __init__(self, graph, index):
        self.graph = graph
        self.index = index
        self.adjacency_list = [] # List of tuples (arc, vertex)
    
    def __str__(self):
        return u'<Vertex %d>' % self.index
    
    def __unicode__(self):
        return self.__str__()


class Arc:
    def __init__(self, graph, vertices, m=3):
        self.graph = graph
        self.vertices = vertices
        self.m = m
    
    def __str__(self):
        return self.vertices.__str__()



class CoxeterGraph:
    SPHERICAL = 'S'
    AFFINE = 'A'
    
    def __init__(self):
        pass
    
    
    def create_arc(self, vertices, m=3):
        arc = Arc(self, vertices, m)
        self.arcs.append(arc)
        
        for i in xrange(2):
            v = vertices[i]
            w = vertices[1-i]
            v.adjacency_list.append((arc, w))
    
    
    def create_graph_from_arcs(self, arcs):
        """
        Fills self.vertices and self.arcs from the list of arcs in the form (i,j,m).
        """
        self.vertices = { v: Vertex(self, v) for v in set([i for (i,j,m) in arcs] + [j for (i,j,m) in arcs]) }
        self.arcs = []
        for (i,j,m) in arcs:
            self.create_arc((self.vertices[i], self.vertices[j]), m)
        
        if self.size == 1:
            # if there is only one vertex, it was not added because there are no arcs
            self.vertices = { 1: Vertex(self, 1) }
    
    
    def find_connected_components(self, simplex):
        """
        Returns the list of connected components of the subgraph induced by the simplex.
        Each connected component is given as a simplex.
        """
        visited = set() # indices of visited vertices
        queue = []
        components = []
        component = []
        
        for i in self.vertices.iterkeys():
            if not i in visited:
                queue.append(i)
                
                while len(queue)>0:
                    j = queue.pop(0)
                    if not j in visited and j in simplex:
                        visited.add(j)
                        component.append(j)
                        for (a,k) in self.vertices[j].adjacency_list:
                            queue.append(k.index)
                
                if len(component) > 0:
                    # found a new component
                    components.append(tuple(sorted(component)))
                    component = []
        
        return components
    
    
    def get_coxeter_type(self, simplex):
        """
        Returns the Coxeter type of the subgraph induced by the simplex (assumed to be connected!).
        This method must be overwritten by the subclass.
        """
        raise Exception()
    
    
    def weight(self, simplex):
        """
        Returns the weight of the simplex.
        """
        components = self.find_connected_components(simplex)
        return sum([Weight()] + [self.get_coxeter_type(component).weight() for component in components])
    
    
    def is_simplex_relevant(self, simplex):
        """
        Says if the simplex is relevant (when deciding if the matching is precise).
        For instance, in A_n with f>0 or g>0, only simplices containing the first f vertices and the last g vertices are relevant.
        This method must be overwritten by the subclass in order to change the default behaviour.
        """
        return True
    
    
    def __str__(self):
        return [map(lambda x: x.index, a.vertices) for a in self.arcs].__str__()




class SphericalACoxeterGraph(CoxeterGraph):
    def __init__(self, n, f=0, g=0):
        self.category = self.SPHERICAL
        self.type = 'A'
        self.n = n
        self.size = n
        # vertices are numbered 1,2,...,n
        assert n >= 0
        
        # create graph
        arcs = [(i, i+1, 3) for i in xrange(1, n)]
        self.create_graph_from_arcs(arcs)
        
        assert 0 <= f <= n
        assert 0 <= g <= n
        self.f = f  # the first f vertices are special
        self.g = g  # the last g vertices are special
    
    def get_coxeter_type(self, simplex):
        """
        Returns the Coxeter type of the subgraph induced by the simplex (assumed to be connected!).
        """
        # all irreducible subgraphs are of type A_n.
        return CoxeterType('A', len(simplex))
    
    def is_simplex_relevant(self, simplex):
        """
        Says if the simplex is relevant (when deciding if the matching is precise).
        Only simplices containing the first f vertices and the last g vertices are relevant.
        """
        return all(i in simplex for i in range(1, self.f+1)+range(self.n-self.g+1, self.n+1))


class SphericalBCoxeterGraph(CoxeterGraph):
    def __init__(self, n, g=0):
        self.category = self.SPHERICAL
        self.type = 'B'
        self.n = n
        self.size = n
        # vertices are numbered 1,2,...,n
        # the special edge with m=4 is (1,2)
        assert n >= 2
        
        # create graph
        arcs = [(1,2,4)] + [(i, i+1, 3) for i in xrange(2, n)]
        self.create_graph_from_arcs(arcs)
        
        assert 0 <= g <= n
        self.g = g  # the last g vertices are special
    
    def get_coxeter_type(self, simplex):
        """
        Returns the Coxeter type of the subgraph induced by the simplex (assumed to be connected!).
        """
        # the subgraph is of type B_n if 1 and 2 are present, A_n otherwise.
        if all(i in simplex for i in [1,2]):
            return CoxeterType('B', len(simplex))
        else:
            return CoxeterType('A', len(simplex))
    
    def is_simplex_relevant(self, simplex):
        """
        Says if the simplex is relevant (when deciding if the matching is precise).
        Only simplices containing the last g vertices are relevant.
        """
        return all(i in simplex for i in xrange(self.n-self.g+1, self.n+1))


class SphericalDCoxeterGraph(CoxeterGraph):
    def __init__(self, n, g=0):
        self.category = self.SPHERICAL
        self.type = 'D'
        self.n = n
        self.size = n
        # vertices are numbered 1,2,...,n
        # the edges are (1,3), (2,3), (3,4), (4,5), ..., (n-1, n)
        assert n >= 4
        
        # create graph
        arcs = [(1,3,3), (2,3,3)] + [(i, i+1, 3) for i in xrange(3, n)]
        self.create_graph_from_arcs(arcs)
        
        assert 0 <= g <= n
        self.g = g  # the last g vertices are special
    
    def get_coxeter_type(self, simplex):
        """
        Returns the Coxeter type of the subgraph induced by the simplex (assumed to be connected!).
        """
        # the subgraph is of type D_n if 1,2,3,4 are present, A_n otherwise.
        if all(i in simplex for i in [1,2,3,4]):
            return CoxeterType('D', len(simplex))
        else:
            return CoxeterType('A', len(simplex))
    
    def is_simplex_relevant(self, simplex):
        """
        Says if the simplex is relevant (when deciding if the matching is precise).
        Only simplices containing the last g vertices are relevant.
        """
        return all(i in simplex for i in xrange(self.n-self.g+1, self.n+1))


class AffineACoxeterGraph(CoxeterGraph):
    def __init__(self, n):
        self.category = self.AFFINE
        self.type = 'A'
        self.n = n
        self.size = n+1
        
        # vertices are numbered 0,1,...,n
        
        # create graph
        arcs = [(i, (i+1)%(n+1), 3) for i in xrange(n+1)]
        self.create_graph_from_arcs(arcs)
    
    def get_coxeter_type(self, simplex):
        """
        Returns the Coxeter type of the subgraph induced by the simplex (assumed to be connected!).
        """
        # all irreducible spherical subgraphs are of type A_n.
        return CoxeterType('A', len(simplex))


class AffineBCoxeterGraph(CoxeterGraph):
    def __init__(self, n):
        self.category = self.AFFINE
        self.type = 'B'
        self.n = n
        self.size = n+1
        
        # vertices are numbered 0,1,...,n
        # the edges are (0,1) with m=4, (1,2), ..., (n-2,n-1), (n-2,n)
        assert n >= 3
        
        # create graph
        arcs = [(0, 1, 4)] + [(i, i+1, 3) for i in xrange(1, n-1)] + [(n-2, n, 3)]
        self.create_graph_from_arcs(arcs)
    
    def get_coxeter_type(self, simplex):
        """
        Returns the Coxeter type of the subgraph induced by the simplex (assumed to be connected!).
        """
        if all(i in simplex for i in [0,1]):
            return CoxeterType('B', len(simplex))
        elif all(i in simplex for i in [self.n, self.n-1, self.n-2, self.n-3]):
            return CoxeterType('D', len(simplex))
        else:
            # all the other irreducible spherical subgraphs are of type A_n.
            return CoxeterType('A', len(simplex))


class AffineCCoxeterGraph(CoxeterGraph):
    def __init__(self, n):
        self.category = self.AFFINE
        self.type = 'C'
        self.n = n
        self.size = n+1
        
        # vertices are numbered 0,1,...,n
        # the edges are (0,1) with m=4, (1,2), ..., (n-2,n-1), (n-1,n) with m=4
        assert n >= 2
        
        # create graph
        arcs = [(0, 1, 4)] + [(i, i+1, 3) for i in xrange(1,n-1)] + [(n-1, n, 4)]
        self.create_graph_from_arcs(arcs)
    
    def get_coxeter_type(self, simplex):
        """
        Returns the Coxeter type of the subgraph induced by the simplex (assumed to be connected!).
        """
        if all(i in simplex for i in [0,1]) or all(i in simplex for i in [self.n-1, self.n]):
            return CoxeterType('B', len(simplex))
        else:
            # all the other irreducible spherical subgraphs are of type A_n.
            return CoxeterType('A', len(simplex))


class AffineDCoxeterGraph(CoxeterGraph):
    def __init__(self, n):
        self.category = self.AFFINE
        self.type = 'D'
        self.n = n
        self.size = n+1
        
        # vertices are numbered 0,1,...,n
        # the edges are (0,2), (1,2), ..., (n-2,n-1), (n-2,n)
        assert n >= 4
        
        # create graph
        arcs = [(0, 2, 3)] + [(i, i+1, 3) for i in xrange(1,n-1)] + [(n-2, n, 3)]
        self.create_graph_from_arcs(arcs)
    
    def get_coxeter_type(self, simplex):
        """
        Returns the Coxeter type of the subgraph induced by the simplex (assumed to be connected!).
        """
        if all(i in simplex for i in [0,1,2,3]) or all(i in simplex for i in [self.n-3,self.n-2,self.n-1,self.n]) or (self.n==4 and len(simplex)==4):
            return CoxeterType('D', len(simplex))
        else:
            return CoxeterType('A', len(simplex))


class SphericalExceptionalCoxeterGraph(CoxeterGraph):
    def __init__(self, type, n, m=None):
        self.category = self.SPHERICAL
        self.type = type
        self.n = n
        self.size = n
        self.m = m
        # vertices are numbered 1,2,...,n
        assert n >= 0
        
        # create graph
        if type == 'E':
            assert 6 <= n <= 8
            arcs = [(i, i+1, 3) for i in xrange(1, n-1)] + [(n-3, n, 3)]
        elif type == 'F':
            assert n == 4
            arcs = [(1, 2, 3), (2, 3, 4), (3, 4, 3)]
        elif type == 'H':
            assert 3 <= n <= 4
            arcs = [(1, 2, 5)] + [(i, i+1, 3) for i in xrange(2, n)]
        elif type == 'I':
            assert n == 2
            arcs = [(1, 2, m)]
        else:
            raise Exception("Invalid parameters for exceptional Coxeter graph")
        
        self.create_graph_from_arcs(arcs)
    
    def get_coxeter_type(self, simplex):
        """
        Returns the Coxeter type of the subgraph induced by the simplex (assumed to be connected!).
        """
        type = self.type
        n = self.n
        
        if type == 'E':
            if all(v in simplex for v in [n, n-1, n-5]):
                # type E
                return CoxeterType('E', len(simplex))
            elif all(v in simplex for v in [n, n-2, n-4]):
                # type D
                return CoxeterType('D', len(simplex))
            else:
                # type A
                return CoxeterType('A', len(simplex))
        
        elif type == 'F':
            if len(simplex) == 4:
                # type F
                return CoxeterType('F', 4)
            elif all(v in simplex for v in [2, 3]):
                # type B
                return CoxeterType('B', len(simplex))
            else:
                # type A
                return CoxeterType('A', len(simplex))
        
        elif type == 'H':
            if all(v in simplex for v in [1, 3]):
                # type H
                return CoxeterType('H', len(simplex))
            elif all(v in simplex for v in [1, 2]):
                # type I_2(5)
                return CoxeterType('I', 2, 5)
            else:
                # type A
                return CoxeterType('A', len(simplex))

        elif type == 'I':
            if len(simplex) == 2:
                # type I_2(m)
                return CoxeterType('I', 2, self.m)
            else:
                # type A
                return CoxeterType('A', len(simplex))



class AffineExceptionalCoxeterGraph(CoxeterGraph):
    def __init__(self, type, n, m=None):
        self.category = self.AFFINE
        self.type = type
        self.n = n
        self.size = n+1
        self.m = m
        # vertices are numbered 0,1,2,...,n
        assert n >= 0
        
        # create graph
        if type == 'tE':
            assert 6 <= n <= 8
            arcs = [(i, i+1, 3) for i in xrange(n-2)]
            if n == 6:
                arcs += [(2, 5, 3), (5, 6, 3)]
            elif n == 7:
                arcs += [(5, 6, 3), (3, 7, 3)]
            elif n == 8:
                arcs += [(6, 7, 3), (5, 8, 3)]
            else:
                raise Exception
        elif type == 'tF':
            assert n == 4
            arcs = [(0, 1, 3), (1, 2, 4), (2, 3, 3), (3, 4, 3)]
        elif type == 'tG':
            assert n == 2
            arcs = [(0, 1, 6), (1, 2, 3)]
        elif type == 'tI':
            assert n == 1
            arcs = [(0, 1, None)]
        else:
            raise Exception("Invalid parameters for exceptional Coxeter graph")
        
        self.create_graph_from_arcs(arcs)
    
    def get_coxeter_type(self, simplex):
        """
        Returns the Coxeter type of the subgraph induced by the simplex (assumed to be connected!).
        """
        type = self.type
        n = self.n
        
        if type == 'tE' and n == 6:
            if len([v for v in simplex if v in [0, 4, 6]]) == 2 and len(simplex) == 6:
                # type E_6
                return CoxeterType('E', 6)
            elif all(v in simplex for v in [1, 3, 5]):
                # type D
                return CoxeterType('D', len(simplex))
            else:
                # type A
                return CoxeterType('A', len(simplex))
        
        elif type == 'tE' and n == 7:
            if all(v in simplex for v in [1, 5, 7]):
                # type E
                return CoxeterType('E', len(simplex))
            elif all(v in simplex for v in [2, 4, 7]):
                # type D
                return CoxeterType('D', len(simplex))
            else:
                # type A
                return CoxeterType('A', len(simplex))
        
        elif type == 'tE' and n == 8:
            if all(v in simplex for v in [3, 7, 8]):
                # type E
                return CoxeterType('E', len(simplex))
            elif all(v in simplex for v in [4, 6, 8]):
                # type D
                return CoxeterType('D', len(simplex))
            else:
                # type A
                return CoxeterType('A', len(simplex))

        elif type == 'tF':
            if all(v in simplex for v in [0, 1, 2, 3]):
                # type F
                return CoxeterType('F', 4)
            if all(v in simplex for v in [1, 2]):
                # type B
                return CoxeterType('B', len(simplex))
            else:
                # type A
                return CoxeterType('A', len(simplex))
        
        elif type == 'tG':
            if all(v in simplex for v in [0, 1]):
                # type I_2(6)
                return CoxeterType('I', 2, 6)
            else:
                # type A
                return CoxeterType('A', len(simplex))
        
        elif type == 'tI':
            # type A
            return CoxeterType('A', len(simplex))
        
        else:
            raise Exception("Invalid parameters for exceptional Coxeter graph")


