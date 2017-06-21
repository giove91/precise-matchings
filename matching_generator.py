#!/usr/bin/python
# coding=utf8

from coxeter_graph import CoxeterGraph

from itertools import combinations, chain, permutations

import copy


def powerset(x):
    """
    Returns the powerset of the iterable x, as a (generator of a) list of tuples.
    """
    return chain.from_iterable(combinations(x,r) for r in xrange(len(x)+1))



class MatchingGenerator:
    
    def __init__(self, debug=False):
        self.clear_matching()
        self.debug = debug
    
    
    def clear_matching(self):
        self.matching = []
    
    
    def add_to_matching(self, sigma, v):
        """
        Add to matching the arc (sigma \cup {v}, sigma \setminus {v}).
        """
        a = list(sigma)
        if not v in sigma:
            a.append(v)
        b = [i for i in a if i != v]
        
        x = tuple(sorted(a))
        y = tuple(sorted(b))
        
        self.matching.append((x,y))
        
        if self.debug:
            print "add to matching", x, y, "\tfrom", sigma
    
    
    def generate_matching(self, simplicial_complex, d, **kwargs):
        """
        Generates the matching for the given simplicial complex (associated to some Coxeter graph).
        """
        self.clear_matching()
        
        coxeter_graph = simplicial_complex.coxeter_graph
        category = coxeter_graph.category
        type = coxeter_graph.type
        n = coxeter_graph.n
        
        for sigma in simplicial_complex.simplices.iterkeys():
            v = None
            if category == CoxeterGraph.SPHERICAL:
                if type == 'A':
                    v = self.generate_spherical_A_matching(n, d, sigma, **kwargs)
                elif type == 'B':
                    v = self.generate_spherical_B_matching(n, d, sigma, **kwargs)
                elif type == 'D':
                    v = self.generate_spherical_D_matching(n, d, sigma, **kwargs)
                else:
                    raise NotImplementedError
            
            elif category == CoxeterGraph.AFFINE:
                if type == 'A':
                    v = self.generate_affine_A_matching(n, d, sigma, **kwargs)
                elif type == 'B':
                    v = self.generate_affine_B_matching(n, d, sigma, **kwargs)
                elif type == 'C':
                    v = self.generate_affine_C_matching(n, d, sigma, **kwargs)
                elif type == 'D':
                    v = self.generate_affine_D_matching(n, d, sigma, **kwargs)
                else:
                    raise NotImplementedError
            
            if v is not None:
                self.add_to_matching(sigma, v)
        
        matching = self.matching
        self.matching = list(set(matching))
        # print sorted(matching), len(matching)
        # print sorted(self.matching), len(self.matching)
        if self.debug:
            print [x for x in self.matching if matching.count(x) != 2]
        assert len(matching) == 2*len(self.matching)
    
    
    ### Finite type ###
    
    def generate_spherical_A_matching_independent(self, n, d, sigma, f=0):
        """
        Generates the matching for the A_n case, leaving only independent subgraphs as critical simplices.
        Vertices are considered to be labeled 1,2,...,n.
        The matching is restricted to the simplices sigma such that 1,2,...,f are in sigma.
        """
        s = set(sigma)
        
        assert n >= 0
        assert 0 <= f <= n
        
        if any(i not in sigma for i in xrange(1, f+1)):
            # sigma does not contain 1,2,...,f-1, so it must not be considered
            return None
        
        elif f == n:
            # nothing can be done
            return None
        
        elif f >= d:
            # recursively generate the matching on the last n-d vertices
            vertex = self.generate_spherical_A_matching_independent(n-d, d, [i-d for i in sigma if i>=d+1], f=f-d)
            if vertex is not None:
                return vertex + d
            else:
                return None
        
        elif f == d-1:
            # perfect matching adding and removing d, if n >= d
            return d if n >= d else None
        
        assert f <= d-2
        
        # search for d-1 consecutive vertices
        v = None
        first = None
        for i in xrange(1,n+1):
            if first is not None and i-first == d-1:
                v = i
                break
        
            if first is None:
                if i in s:
                    first = i
            else:
                if i not in s:
                    first = None
        
        if v is not None:
            return v
    
    
    def generate_spherical_A_matching(self, n, d, sigma, f=0, g=0):
        """
        Generates the matching for the A_n case, also involving 0-weight simplices.
        Vertices are considered to be labeled 1,2,...,n.
        The matching is restricted to the simplices sigma such that the vertices 1,...,f,n-g+1,...,n are in sigma.
        """
        s = set(sigma)
        
        assert n >= 0
        
        assert 0 <= f <= n and 0 <= g <= n
        
        if any(i not in sigma for i in range(1, f+1)+range(n-g+1,n+1)):
            # sigma does not contain 1,2,...,f, so it must not be considered
            return None
        
        elif f+g >= n:
            # nothing can be done
            return None
        
        
        elif f >= d:
            # recursively generate the matching on the last n-d vertices
            vertex = self.generate_spherical_A_matching(n-d, d, [i-d for i in sigma if i>=d+1], f=f-d, g=g)
            if vertex is not None:
                return vertex + d
            else:
                return None
        
        elif g >= d:
            # recursively generate the matching on the first n-d vertices
            vertex = self.generate_spherical_A_matching(n-d, d, [i for i in sigma if i<=n-d], f=f, g=g-d)
            if vertex is not None:
                return vertex
            else:
                return None
        
        
        elif f == d-1:
            # perfect matching adding and removing d
            return d
        
        elif g == d-1:
            # perfect matching adding and removing n-d+1
            return n-d+1
        
        assert 0 <= f <= d-2
        
        if n < d+g:
            # in this case I cannot continue recursively!
            # try to add/remove f+1
            if all(i in s for i in xrange(f+2, n+1)) and d-1 <= n <= d-1+f:
                # this is the only case that goes wrong (we have two critical cells)
                return None
            else:
                # add/remove f+1
                return f+1
        
        assert n >= d+g
        
        # find the size of the connected component of 1
        k = min(i for i in xrange(1, n+2) if i not in s) - 1
        
        if k >= d-1:
            # long component
            
            if d in sigma and d <= n-g:
                # remove d (the condition d <= n-g is now automatically verified)
                return d
            
            if d not in sigma and d <= n:
                # add d
                return d
        
        else:
            # short component
            
            if f+1 in s:
                # remove f+1
                assert f+1 <= n-g
                return f+1
            
            else:
                if all(i in s for i in xrange(f+2, d)):
                    # I cannot add f+1 because I would created a too big component
                    # an A_{n-f-1} remains, with f=d-2-f
                    vertex = self.generate_spherical_A_matching(n-f-1, d, [i-f-1 for i in sigma if i>f+1], f=d-2-f, g=g)
                    return vertex + f + 1 if vertex is not None else None
                
                else:
                    # add f+1
                    return f+1
    
    
    def generate_spherical_B_matching(self, n, d, sigma, g=0):
        """
        Generates the matching for the B_n case, also involving 0-weight simplices.
        Vertices are considered to be labeled 1,2,...,n.
        The special edge with m=4 is (1,2).
        The matching is restricted to the simplices sigma such that the vertices n-g+1,...,n are in sigma (this works only for d even).
        """
        if n <= 1:
            return self.generate_spherical_A_matching(n, d, sigma, g=g)
        
        s = set(sigma)
        
        if any(i not in sigma for i in xrange(n-g+1,n+1)):
            # sigma does not contain 1,2,...,f, so it must not be considered
            assert d%2 == 0
            return None
        
        
        if d%2 == 1:
            return 1
        
        else:
            # size of the B_k component
            k = min(i-1 for i in xrange(1, n+1) if i not in s) if len(s) < n else n
            
            if d == 2:
                # try to match with a cell with the same k
                if k < n:
                    vertex = self.generate_spherical_A_matching(n-k-1, d, [i-k-1 for i in sigma if i >= k+2], g=g)
                    if vertex is not None:
                        return vertex + k + 1
            
            else:
                assert d >= 4
                
                # euclidean division
                q = k/(d/2)
                r = k - q*(d/2)
                
                v = q*(d/2) + 1
                
                if r >= 1:
                    # remove v, unless we are in the bad case
                    if not (g > 0 and v > n-g):
                        return v
                    
                    else:
                        return None
                
                else:
                    assert r == 0
                    
                    if v > n:
                        # vertex v does not exist
                        return None
                    
                    elif all(i in s for i in xrange(q*(d/2)+2, (q+1)*(d/2)+1)):
                        # there is a big component (length >= d/2-1) immediately after v,
                        # so I can't add v
                        # an A_{n-d/2-1} remains, with f=d/2-1 and g=g
                        vertex = self.generate_spherical_A_matching(n-q*d/2-1, d, [i-q*d/2-1 for i in sigma if i>q*d/2+1], f=d/2-1, g=g)
                        return vertex + q*d/2 + 1 if vertex is not None else None
                    
                    else:
                        # add v
                        return v
    
    
    def generate_spherical_D_matching(self, n, d, sigma, g=0):
        """
        Generates the matching for the D_n case.
        Vertices are considered to be labeled 1,2,...,n.
        The edges are (1,3), (2,3), (3,4), (3,5), ...
        The matching is restricted to the simplices sigma such that the vertices n-g+1,...,n are in sigma (this works only for d even).
        """
        s = set(sigma)
        
        if any(i not in sigma for i in xrange(n-g+1,n+1)):
            # sigma does not contain n-g+1,...,n, so it must not be considered
            assert d%2 == 0
            return None
        
        
        if n <= 1:
            # A_n
            return self.generate_spherical_A_matching(n, d, sigma, g=g)
        elif n == 2:
            # two copies of A_1
            if d == 2 or g >= 2:
                # no matching can be done
                return None
            else:
                # perfect matching adding and removing 1
                return 1
        elif n == 3:
            # A_3
            if g != 1:
                new_vertices = {1: 1, 2: 3, 3: 2}
                new_vertices_inv = {i: j for (j,i) in new_vertices.iteritems()}
                vertex = self.generate_spherical_A_matching(3, d, sorted([new_vertices_inv[i] for i in sigma]), g=g)
                return new_vertices[vertex] if vertex is not None else None
            else:
                assert g == 1
                # for g=1, the fixed vertex of A_3 is the one in the middle!
                # this case should be done separately
                if 2 not in s and d != 3:
                    return 1
                elif 2 in s and d not in [2, 4]:
                    return 1
                else:
                    return None
        
        
        if d%2 == 1:
            # d odd
            if not 1 in s:
                # 1 not in sigma
                vertex = self.generate_spherical_A_matching(n-1, d, [i-1 for i in sigma])
                if vertex is not None:
                    return vertex + 1
            else:
                # perfect matching adding and removing 2
                return 2
        
        else:
            if 2 not in s:
                # this case is in common between d=2 and d>=4
                # a A_{n-1} remains
                new_vertices = {i: n+1-i for i in xrange(1,n-1)}
                new_vertices[n-1] = 1
                new_vertices_inv = {i: j for (j,i) in new_vertices.iteritems()}
                vertex = self.generate_spherical_A_matching(n-1, d, sorted([new_vertices_inv[i] for i in sigma]), f=g)
                return new_vertices[vertex] if vertex is not None else None
            
            elif d == 2 and any(i not in s for i in [1,2,3,4]):
                if all(i in s for i in [1,2,4]):
                    assert 3 not in s
                    if 5 in s and 5 <= n-g:
                        # remove 5
                        return 5
                    elif 5 not in s and n >= 5:
                        # add 5
                        return 5
                    else:
                        return None
                
                else:
                    # 2 in sigma, at least one of 1, 4 not in sigma
                    assert 2 in s
                    assert any(i not in s for i in [1,4])
                    # try to add and remove 3
                    if 3 in s and 3 <= n-g:
                        # remove 3
                        return 3
                    elif 3 not in s:
                        # add 3
                        return 3
                    else:
                        return None
            
            elif d >= 4 and 3 not in s:
                # perfect matching adding and removing 1
                return 1
            
            elif d >= 4 and 4 not in s:
                if d == 4:
                    # two A_{n-4} remain (depending on whether 1 is in sigma or not)
                    vertex = self.generate_spherical_A_matching(n-4, d, sorted([n+1-i for i in s if i>4]), f=g)
                    return n+1 - vertex if vertex is not None else None
                else:
                    # perfect matching adding and removing 1
                    return 1
            
            elif d >= 4 and 1 not in s:
                assert 4 in s
                
                if all(i in s for i in xrange(2, d/2+2)):
                    # a A_{n-1} with g=max(d/2, 3) remains
                    vertex = self.generate_spherical_A_matching(n-1, d, sorted([n+1-i for i in sigma]), f=g, g=max(d/2, 3))
                    return n+1 - vertex if vertex is not None else None
                else:
                    # add 1
                    # (if the component of 2 has at most d/2-1 vertices, then it is matched with one of the simplices below)
                    return 1
            
            else:
                # this case is in common between d=2 and d>=4
                assert all(i in s for i in [1,2,3,4])
                k = min(i for i in xrange(1, n+1) if i not in s) - 1 if len(s) < n else n  # size of the D_k component on the left
                
                # euclidean division
                q = k/(d/2)
                r = k - q*(d/2)
                
                # the case q odd and r=0 should be treated together with q even
                if q%2 == 1 and r==0:
                    q -= 1
                    r = d/2
                
                # find the vertex I want to add/remove, without changing the weight of the D_k component
                v = q*d/2+1 if q%2==0 else q*d/2+2
                
                if v in s:
                    # simply remove v
                    if 2 <= v <= 4:
                        # removing 2,3,4 is dangerous, because there isn't anymore a D_k component!
                        # this actually should not happen!
                        assert False
                        return None
                    elif v > n-g:
                        # cannot remove v
                        return None
                    else:
                        return v
                
                else:
                    if v > n:
                        # the vertex v does not exist
                        return None
                    
                    elif q%2 == 0 and all(i in s for i in xrange(q*d/2+2, (q+1)*d/2+2)):
                        # q is even, and there is a too big component
                        # an A_{n-q*d/2-1} remains, with g=d/2
                        vertex = self.generate_spherical_A_matching(n-q*d/2-1, d, sorted([n+1-i for i in s if i>v]), f=g, g=d/2)
                        return n+1 - vertex if vertex is not None else None
                    
                    elif q%2 == 1 and all(i in s for i in xrange(q*d/2+3, (q+1)*d/2+1)):
                        # q is odd, and there is a too big component
                        # an A_{n-q*d/2-2} remains, with g=d/2-2
                        vertex = self.generate_spherical_A_matching(n-q*d/2-2, d, sorted([n+1-i for i in s if i>v]), f=g, g=d/2-2)
                        return n+1 - vertex if vertex is not None else None
                    
                    else:
                        # add v
                        return v
    
    
    ### Affine type ###
    
    def generate_affine_A_matching(self, n, d, sigma):
        """
        Generates the matching for the tilde A_n case.
        Vertices are considered to be labeled 0,1,...,n.
        """
        s = set(sigma)
        
        # search for the first vertex not in sigma
        k = min([i for i in xrange(n+1) if i not in s])
        
        # if we fix k, an A_n remains with k fixed vertices
        vertex = self.generate_spherical_A_matching(n, d, sorted([(k-i)%(n+1) for i in sigma]), f=k)
        return (k - vertex) % (n+1) if vertex is not None else None
    
    
    
    def generate_affine_B_matching(self, n, d, sigma):
        """
        Generates the matching for the tilde B_n case.
        Vertices are considered to be labeled 0,1,2,...,n.
        """
        s = set(sigma)
        
        if d%2 == 1:
            # almost perfect matching adding and removing 0
            if sigma != tuple(range(1, n+1)):
                return 0
        
        else:
            # size of the B_k component
            k = min(i for i in xrange(n+1) if i not in s)
            
            special_case = False
            if k == n-1 and n in s:
                # sigma=(0,1,...,n-2,n) has k=n (not n-1)!
                k = n
                special_case = True
            
            # Euclidean division
            q = k / (d/2)
            r = k - q*d/2
            
            if r >= 1:
                # idea: remove q*d/2
                if special_case and q*d/2 == n-1:
                    # actually we want to remove n, not n-1, but...
                    # (0,1,...,n-2,n-1) wants to do the same!
                    return None
                return q*d/2
            
            assert r == 0
            if all(i in s for i in xrange(q*d/2+1, (q+1)*d/2)):
                # cannot add k==q*d/2
                vertex = self.generate_spherical_D_matching(n-k, d, sorted([n-i+1 for i in sigma if i>k]), g=d/2-1)
                return n-vertex+1 if vertex is not None else None
            
            else:
                # idea: add k==q*d/2
                if all(i in s for i in xrange(q*d/2+1, (q+1)*d/2-1)) and (q+1)*d/2 in s and (q+1)*d/2 == n:
                    # sigma=(0,1,...,k-1,k+1,...,n-2,n), and adding k would create a (too big) B_n component
                    assert k == n-d/2
                    return None
                
                if len(sigma) == n:
                    # adding k we would obtain (0,1,...,n), impossible
                    return None
            
                return k
    
    
    def generate_affine_C_matching(self, n, d, sigma):
        """
        Generates the matching for the tilde C_n case.
        Vertices are considered to be labeled 0,1,2,...,n.
        """
        s = set(sigma)
        
        if d%2 == 1:
            # almost perfect matching adding and removing 0
            if sigma != tuple(range(1,n+1)):
                return 0
        
        else:
            k = max(j for j in xrange(n+1) if j not in sigma)
            # on the right there is some B_h component, and k vertices remain on the left
            
            if k >= 2:
                # a B_k remains on the left
                vertex = self.generate_spherical_B_matching(k, d, [i+1 for i in sigma if i < k])
                if vertex is not None:
                    return vertex - 1
            else:
                # an A_1 or A_0 remains
                vertex = self.generate_spherical_A_matching(k, d, [i+1 for i in sigma if i < k])
                if vertex is not None:
                    return vertex - 1
    
    
    def generate_affine_D_matching(self, n, d, sigma):
        """
        Generates the matching for the tilde D_n case.
        Vertices are considered to be labeled 0,1,2,...,n.
        """
        s = set(sigma)
        
        if d%2 == 1:
            if all(i in s for i in [1,2,3]):
                # 1,2,3 in sigma
                if not all(i in s for i in xrange(1,n+1)):
                    # perfect matching except for the critical cell (1,2,...,n)
                    return 0
            
            elif 1 not in s:
                # 1 not in sigma
                # a D_n remains
                new_vertices = {i: n+1-i for i in xrange(1,n)}
                new_vertices[n] = 0
                new_vertices_inv = {i: j for (j,i) in new_vertices.iteritems()}
                vertex = self.generate_spherical_D_matching(n, d, [new_vertices_inv[i] for i in sigma])
                if vertex is not None:
                    return new_vertices[vertex]
            
            else:
                # perfect matching adding and removing 0
                return 0
        
        else:
            if n == 4 and d <= 6:
                # n=4 is a degenerate case
                
                if d == 2:
                    if len(s) == 1 and 2 not in s:
                        # add 2
                        return 2
                    elif len(s) == 2:
                        # add/remove 2
                        return 2
                    elif len(s) == 3 and 2 in s:
                        # remove 2
                        return 2
                    else:
                        return None
                
                elif d == 4:
                    if 2 not in s or all(i not in s for i in [1,3,4]):
                        # add/remove 0
                        return 0
                    else:
                        return None
                
                elif d == 6:
                    if 2 in s and 0 not in s and len(s)>=3:
                        return None
                    elif 2 in s and 0 in s and len(s)==4:
                        return None
                    else:
                        # add/remove 0
                        return 0
            
            elif 1 not in s:
                # this case is in common between d=2 and d>=4
                # a D_n remains
                new_vertices = {i: n+1-i for i in xrange(1,n)}
                new_vertices[n] = 0
                new_vertices_inv = {i: j for (j,i) in new_vertices.iteritems()}
                vertex = self.generate_spherical_D_matching(n, d, sorted([new_vertices_inv[i] for i in sigma]))
                return new_vertices[vertex] if vertex is not None else None
            
            elif d == 2 and any(i not in s for i in [0,1,2,3]):
                if all(i in s for i in [0,1,3]):
                    assert 2 not in s
                    # try to add and remove 4
                    if not all(i in s for i in xrange(5,n+1)):
                        return 4
                    else:
                        return None
                
                else:
                    # 1 in sigma, at least one of 0, 3 not in sigma
                    assert 1 in s
                    assert any(i not in s for i in [0,3])
                    # try to add and remove 2
                    if not all(i in s for i in xrange(1,n+1) if i != 2):
                        return 2
                    else:
                        return None
            
            elif d >= 4 and 0 not in s:
                if all(i in s for i in xrange(1, d/2+1)):
                    # a D_n with g=d/2 remains
                    vertex = self.generate_spherical_D_matching(n, d, sorted([n+1-i for i in sigma]), g=d/2)
                    return n+1 - vertex if vertex is not None else None
                elif d/2+1 == n and all(i in s for i in range(1,d/2)+[n]):
                    # the component goes down: sigma = (1,2,3,...,n-2,n), and the component is too long
                    assert n-1 not in s # otherwise the previous "if" would have been true
                    # the simplex (1,2,3,...,n-2,n-1,n) prefers (1,2,3,...,n-2,n-1), so there is no chance for us
                    return None
                elif len(sigma) < n:
                    # add 0, except when sigma=(1,2,...,n)
                    # (if the component of 1 has at most d/2-1 vertices, then it is matched with one of the simplices below)
                    return 0
                else:
                    return None
            
            elif d >= 4 and 2 not in s:
                # perfect matching adding and removing 0
                return 0
            
            elif d >= 4 and 3 not in s:
                if d == 4:
                    # a D_{n-3} remains
                    vertex = self.generate_spherical_D_matching(n-3, d, sorted([n+1-i for i in s if i>3]))
                    return n+1 - vertex if vertex is not None else None
                else:
                    # perfect matching adding and removing 0
                    return 0
            
            
            else:
                # this case is in common between d=2 and d>=4
                assert all(i in s for i in [0,1,2,3])
                k = min(i for i in xrange(n+1) if i not in s)   # size of the D_k component on the left
                if k == n-1 and n in s:
                    # sigma = (0,1,...,n-2,n)
                    k = n
                
                # euclidean division
                q = k/(d/2)
                r = k - q*(d/2)
                
                # the case q odd and r=0 should be treated together with q even
                if q%2 == 1 and r == 0:
                    q -= 1
                    r = d/2
                
                # find the vertex I want to add/remove, without changing the weight of the D_k component
                v = q*d/2 if q%2==0 else q*d/2+1
                
                if d == 4 and q%2 == 1 and r == 1:
                    # there is only one k with this weight, so no such vertex v can exist
                    # simply do recursively what remains
                    if k <= n-2:
                        # a D_{n-k} remains
                        vertex = self.generate_spherical_D_matching(n-k, d, sorted([n+1-i for i in s if i>k]))
                        return n+1 - vertex if vertex is not None else None
                    else:
                        # there isn't enough space to do anything
                        return None
                
                elif k == n and n-1 not in s and v >= n-1:
                    # sigma = (0,1,...,n-2,n), v = n-1 or v = n (should be n or n+1, respectively)
                    # cannot add v+1 (in the first case because (0,1,...,n-2) is already matched with (0,1,...,n-2,n-1),
                    # in the second case because n+1 does not exist)
                    return None
                
                elif v in s:
                    # simply remove v
                    if 1 <= v <= 3:
                        # removing 1,2,3 is dangerous, because there isn't anymore a D_k component!
                        # this should not happen
                        assert False
                        return None
                    else:
                        return v
                
                else:
                    assert v not in s
                    
                    if v > n:
                        # the vertex v does not exist
                        return None
                    
                    # find the connected component c right after v
                    c_max = min(i for i in xrange(v+1, n+2) if i not in s)  # first vertex not in c
                    c_goes_up = False   # n-1 in c
                    c_goes_down = False # n in c
                    if c_max == n+1:
                        # c contains n-1 and n
                        c_goes_up = True
                        c_goes_down = True
                    elif c_max == n:
                        # c contains n-1 but not n
                        c_goes_up = True
                        c_goes_down = False
                    elif c_max == n-1 and n in s:
                        # c contains n but not n-1
                        c_goes_up = False
                        c_goes_down = True
                        c_max += 1  # change c_max so that c_length is computed correctly
                    c_length = c_max - v - 1    # length of c
                    
                    # find the minimum size of c for which we cannot add v
                    limit_size = d/2 if q%2 == 0 else d/2-2
                    
                    if c_length < limit_size and not all(x for x in [c_goes_up, c_goes_down]):
                        # if the length of c is small enough, and c doesn't contain n-1 and n, we can add v
                        return v
                    
                    elif c_goes_down and not c_goes_up and c_length == limit_size:
                        # this case is degenerate, and it shouldn't be matched with anything
                        return None
                    
                    else:
                        # a D_{n-v} remains, with g=limit_size
                        vertex = self.generate_spherical_D_matching(n-v, d, sorted([n+1-i for i in s if i>v]), g=limit_size)
                        return n+1 - vertex if vertex is not None else None


