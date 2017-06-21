from nzmath import matrix, vector # http://tnt.math.se.tmu.ac.jp/nzmath/

class Cell:
    """
    A cell of a complex.
    """
    
    def __init__(self, d, id=-1, label=None):
        self.d = d # Dimension
        self.id = id
        self.label = label # Additional label for the cell
        
        # The following lists are filled in Complex.__init__()
        self.subcells = [] # Adjacency list of sub-cells (as Edges)
        self.supercells = [] # Adjacency list of super-cells (as Edges)
        
        # The following properties are used in Discrete Morse Theory reduction
        self.matching_edge = None
        self.visited = False # Used in DAG traversal
        self.closed = False # Used in DAG traversal
        self.twin = None # The corresponding cell in the new (reduced) complex
        self.aggregate_weight = None # Used to compute the degree of the new edges
    
    def is_matched(self):
        return self.matching_edge is not None
    
    def children(self):
        """
        Return an iterable of tupes (cell, degree), children in the DAG (taking into account the Discrete Morse Theory matching)
        """
        for e in self.subcells:
            if e != self.matching_edge:
                yield (e.low, e.deg)
        if self.matching_edge is not None and self.matching_edge.low == self:
            yield (self.matching_edge.high, self.matching_edge.deg)
    
    def restricted_children(self, k):
        """
        As self.children(), but it returns only edges between k-dimensional cells and (k-1)-dimensional cells
        """
        if self.d == k:
            # Edges going down
            for e in self.subcells:
                if e != self.matching_edge:
                    yield (e.low, e.deg)
        
        if self.d == k-1:
            # Edge going up (if exists)
            if self.is_matched():
                e = self.matching_edge
                if e.high.d == k:
                    yield (e.high, e.deg)
    
    def __repr__(self):
        if self.label is None:
            return '<Cell (%d,%d)>' % (self.d, self.id)
        else:
            return '<Cell (%d,%d) label ' % (self.d, self.id) + self.label.__str__() + '>'
    
    def __hash__(self):
        return hash((self.d, self.id))


class Edge:
    """
    An edge of a complex.
    """
    
    def __init__(self, high, low, deg, matchable=True):
        self.high = high # The higher-dimensional cell
        self.low = low # The lower-dimensional cell
        self.deg = deg # The incidence degree (assumed to be an integer)
        self.matchable = matchable  # Can be used to prevent adding to matching
        
        # The following properties are used in Discrete Morse Theory reduction
        self.is_in_matching = False
    
    def is_matchable(self):
        # Check if it is possible to add this edge to the matching
        if not self.matchable:
            return False
        for c in [self.high, self.low]:
            if c.is_matched():
                return False
        assert not self.is_in_matching # This condition should be automatically verified at this stage
        if not self.deg in [-1, 1]:
            # Only regular edges can be collapsed
            return False
        return True
    
    def add_to_matching(self):
        # Add this edge to the matching
        assert self.is_matchable()
        for c in [self.high, self.low]:
            c.matching_edge = self
        self.is_in_matching = True
    
    def remove_from_matching(self):
        # Remove this edge from the matching
        assert self.is_in_matching
        for c in [self.high, self.low]:
            c.matching_edge = None
        self.is_in_matching = False
    
    def __repr__(self):
        if self.high.label is None or self.low.label is None:
            return '<Edge (%d,%d)->(%d,%d) of degree %d>' % (self.high.d, self.high.id, self.low.d, self.low.id, self.deg)
        else:
            return '<Edge ' + self.high.label.__str__() + '->' + self.low.label.__str__() + ' of degree %d>' % self.deg


class Complex:
    """
    A complex. Only the incidence degree between cells are stored.
    """
    
    def __init__(self, cells, edges, cells_by_dimension=False):
        if len(cells) == 0:
            self.d = -1
            self.cells = {}
        
        elif cells_by_dimension:
            # The cells are already given by dimension
            self.d = max(cells.iterkeys())
            self.cells = cells
            while self.d in self.cells and len(self.cells[self.d]) == 0:
                # There are no cells of dimension self.d
                del self.cells[self.d]
                self.d -= 1
            if len(self.cells) == 0:
                self.d = -1
        
        else:
            self.d = max(c.d for c in cells) # Dimension
            
            self.cells = {k: [] for k in xrange(min(c.d for c in cells), self.d+1)}
            for c in cells:
                self.cells[c.d].append(c)
        
        # Set the id of the cells
        for l in self.cells.itervalues():
            for (i,c) in enumerate(l):
                c.id = i
        
        self.edges = edges
        
        # Add edges to the adjacency lists
        for e in edges:
            e.high.subcells.append(e)
            e.low.supercells.append(e)
    
    def all_cells(self):
        # Returns an iterable with all the cells (not by dimension)
        for l in self.cells.itervalues():
            for c in l:
                yield c
    
    
    def DFS_acyclic(self, cell, k, print_cycle=False):
        """
        DFS for the purpose of testing acyclicity
        """
        if cell.visited:
            if cell.closed:
                return True
            else:
                if print_cycle:
                    print cell
                return False
        
        # print "Visiting cell (%d,%d)" % (cell.d, cell.id)
        cell.visited = True
        
        # print "Restricted children (%d):" % k, [c for c in cell.restricted_children(k)]
        for (c,w) in cell.restricted_children(k):
            # print "Launching visit of", c, "from", cell
            res = self.DFS_acyclic(c, k, print_cycle)
            if not res:
                if print_cycle:
                    print cell
                return False
        
        cell.closed = True
        return True
    
    def is_acyclic(self, starting_cell, k, print_cycle=False):
        """
        Test acyclicity for k-dimensional and (k-1)-dimensional cells.
        """
        # Initialize DFS
        for c in self.cells[k] + self.cells[k-1]:
            c.visited = False
            c.closed = False
        
        # Visit
        return self.DFS_acyclic(starting_cell, k, print_cycle)
    
    
    def DFS_weight(self, cell, target, k):
        """
        DFS for the purpose of finding weights of the new edges (via dynamic programming on the DAG)
        """
        if cell == target:
            return [1,0]
        
        if cell.visited:
            return cell.aggregate_weight
        
        cell.visited = True
        cell.aggregate_weight = [0,0] # (weight of path from cell to target, passing by an even number of other k-dimensional cells; odd)
        
        for (c,w) in cell.restricted_children(k):
            (x,y) = self.DFS_weight(c, target, k)
            if cell.d == k:
                cell.aggregate_weight[0] += w*x
                cell.aggregate_weight[1] += w*y
            else:
                cell.aggregate_weight[0] += w*y
                cell.aggregate_weight[1] += w*x
        
        return cell.aggregate_weight
    
    
    def morse_reduction(self):
        """
        Perform Discrete Morse Theory collapses, returning a smaller complex with the same homotopy type.
        """
        # Create new cells
        new_cells = {k: [] for k in self.cells.iterkeys()}
        for c in self.all_cells():
            if not c.is_matched():
                c1 = Cell(d=c.d, label=c.label)
                new_cells[c.d].append(c1)
                c1.twin = c
                c.twin = c1
        
        # Create new edges
        new_edges = []
        for k in sorted(self.cells.iterkeys())[1:]:
            # Create edges from k-dimensional cells to (k-1)-dimensional cells
            for a in new_cells[k]:
                for b in new_cells[k-1]:
                    # Initialize DFS
                    for c in self.cells[k] + self.cells[k-1]:
                        c.visited = False
                    
                    # Visit
                    weight = self.DFS_weight(a.twin, b.twin, k)
                    w = weight[0] - weight[1]
                    if w != 0:
                        new_edges.append(Edge(a, b, w))
        
        # Return the new complex
        return Complex(new_cells, new_edges, cells_by_dimension=True)
    
    
    def get_boundaries(self):
        """
        Compute boundary matrices.
        """
        ranks = {k: len(l) for (k, l) in self.cells.iteritems()}
        boundaries = {}
        
        for k in ranks.iterkeys():
            if k-1 not in ranks:
                continue
            
            # Initialize 0-column matrix
            numrows = max(ranks[k-1], 1) # The max with 1 is required by NZMATH
            delta = matrix.Matrix(row=numrows, column=1)
            delta.deleteColumn(1)
            
            for c in self.cells[k]:
                column = vector.Vector([0]*numrows)
                for e in c.subcells:
                    # Update the column of the boundary map
                    c1 = e.low
                    deg = e.deg
                    column[c1.id+1] = deg
                
                # Insert column in the boundary map
                delta.extendColumn(column)
            
            # Adjust the case ranks[k-1]==0
            if ranks[k-1] == 0:
                assert delta.row == 1
                delta.deleteRow(1)
            
            boundaries[k] = delta
        
        return ranks, boundaries


