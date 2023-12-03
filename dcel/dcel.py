#!/usr/bin/env python
# Copyright 2008, Angel Yanguas-Gil

__all__ = ['Dcel', 'Vertex', 'Hedge', 'Face']

from dcel.xygraph import Xygraph
import numpy as np
import math as m


class DcelError(Exception): pass


class Vertex:
    """Minimal implementation of a vertex of a 2D dcel"""

    def __init__(self, px, py):
        """
        initialize the attributes of a vertex.
        :param px: x coordinate of p.
        :param py: y coordinate of p.
        :attr hedgelist: list of half edges originating at p.
        """
        self.x = px
        self.y = py
        self.hedgelist = []

    def sortincident(self):
        """
        sort the out half-edges in anti-clockwise by incident point.
        :return: NA
        """
        self.hedgelist.sort(key=(lambda b: b.angle), reverse=True)


class Hedge:
    """Minimal implementation of a half-edge of a 2D dcel"""

    def __init__(self, v1, v2):
        """
        initializes the vertices of a half-edge.
        :param v1: v1 is the incident vertex of the type Vertex.
        :param v2: v2 is the origin vertex.
        :attr origin: v2.
        :attr twin: the half-edge going from v1 to v2
        :attr face: the associated face of the half-edge
        :attr nexthedge: the next half-edge at the boundary of its face
        :attr prevhedge: the previous half edge at the boundary of its face
        :attr length: the length of the half-edge
        :attr legal: whether the half-edge is legal or not, in context of lawson edge swap.
        :attr external: whether the half-edge is at the boundary of the face at infinity.

        """
        # The origin is defined as the vertex it points to
        self.origin = v2
        self.twin = None
        self.face = None
        self.nexthedge = None
        self.angle = hangle(v2.x - v1.x, v2.y - v1.y)
        self.prevhedge = None
        self.length = m.sqrt((v2.x - v1.x) ** 2 + (v2.y - v1.y) ** 2)
        self.legal = False
        self.external = None


class Face:
    """Implements a face of a 2D dcel"""

    def __init__(self):
        """
        initializes the face
        :attr wedge: a half-edge in the boundary of the face
        :attr data: data corresponding to the face
        :attr external: whehter the face is the face at infinity
        """
        self.wedge = None
        self.data = None
        self.external = None

    def area(self):  # using the shoelace formula for finding area of an arbitrary polygon with vertices given
        """
        find the area of the face using the shoe-lace formula.
        :return:
            area of the face.
        """
        h = self.wedge
        a = 0
        while (h.nexthedge != self.wedge):
            p1 = h.origin
            p2 = h.nexthedge.origin
            a += p1.x * p2.y - p2.x * p1.y
            h = h.nexthedge

        p1 = h.origin
        p2 = self.wedge.origin
        a = (a + p1.x * p2.y - p2.x * p1.y) / 2
        return a

    def perimeter(self):
        """
        find the perimeter of the face by iterating over its edges.
        :return:
            perimeter of the face.
        """
        h = self.wedge
        p = 0
        while (h.nexthedge != self.wedge):
            p += h.length
            h = h.nexthedge
        p += h.length
        return p

    def vertexlist(self):
        """
        find the list of vertices making the face
        :return:
            return the list of vertices making the face.
        """
        h = self.wedge
        pl = [h.origin]
        while (h.nexthedge != self.wedge):
            h = h.nexthedge
            pl.append(h.origin)
        return pl

    def isinside(self, p):
        """
        Determines whether a point is inside a face
        :return:
            detects is a point is inside a face or not.
        """

        h = self.wedge
        inside = False
        if lefton(h, p):
            while (h.nexthedge != self.wedge):
                h = h.nexthedge
                if not lefton(h, p):
                    return False
            return True
        else:
            return False


class Dcel(Xygraph):
    """
    Implements a doubly-connected edge list
    """

    def __init__(self, vl=[], el=[]):
        """
        initializes the dcel structure.
        :param vl: list of vertices as a list of lists
        :param el: list of edges as a list of lists indicating the indices of the points being joined from the vertex list.

        :attr vertices: list of vertex objects in the dcel
        :attr hedges: list of hedges objects in the dcel

        :attr faces: list of face objects in the dcel :attr upper_chain: the list of vertices in x-increasing order that
        form the upper part of the polygon when cut at x-extreme points

        :attr lower_chain: the list of vertices in x-increasing order that form the lower
        part of the polygon when cut at x-extreme points
        """
        Xygraph.__init__(self, vl, el)
        self.vertices = []
        self.hedges = []
        self.faces = []
        self.upper_chain = None
        self.lower_chain = None
        if vl != []:
            self.build_dcel()

    def build_dcel(self):
        """
        Creates the dcel from the list of vertices and edges
        """

        # Step 1: vertex list creation
        for v in self.vl:
            self.vertices.append(Vertex(v[0], v[1]))

        # Step 2: hedge list creation. Assignment of twins and
        # vertices
        for e in self.el:
            if e[0] >= 0 and e[1] >= 0:
                h1 = Hedge(self.vertices[e[0]], self.vertices[e[1]])  # origin is self.vertices[e[1]]
                h2 = Hedge(self.vertices[e[1]], self.vertices[e[0]])  # origin is self.vertices[e[0]]
                h1.twin = h2
                h2.twin = h1
                self.vertices[e[1]].hedgelist.append(h1)  # out going edge
                self.vertices[e[0]].hedgelist.append(h2)  # out going edge
                self.hedges.append(h2)
                self.hedges.append(h1)

        # Step 3: Identification of next and prev hedges
        for v in self.vertices:
            v.sortincident()
            l = len(v.hedgelist)
            if l < 2:
                raise DcelError(
                    "Badly formed dcel: less than two hedges in vertex")
            else:
                for i in range(l - 1):
                    v.hedgelist[i].twin.nexthedge = v.hedgelist[i + 1]
                    v.hedgelist[i + 1].prevhedge = v.hedgelist[i].twin
                v.hedgelist[l - 1].twin.nexthedge = v.hedgelist[0]
                v.hedgelist[0].prevhedge = v.hedgelist[l - 1].twin

        # Step 4: Face assignment
        provlist = self.hedges[:]
        nf = 0
        nh = len(self.hedges)

        while nh > 0:
            h = provlist.pop()
            nh -= 1
            # We check if the hedge already points to a face
            if h.face == None:
                f = Face()
                nf += 1
                # We link the hedge to the new face
                f.wedge = h
                f.wedge.face = f
                # And we traverse the boundary of the new face
                while (h.nexthedge != f.wedge):
                    h = h.nexthedge
                    h.face = f
                self.faces.append(f)
        # And finally we have to determine the external face
        for f in self.faces:
            f.external = f.area() < 0

        for h in self.hedges:
            h.external = h.face.external
            if h.face.external:
                h.legal = True
                h.twin.legal = True

    def findpoints(self, pl, onetoone=False):
        """Given a list of points pl, returns a list of
        with the corresponding face each point belongs to and
        None if it is outside the map.

        """

        ans = []
        if onetoone:
            fl = self.faces[:]
            for p in pl:
                found = False
                for f in fl:
                    if f.external:
                        continue
                    if f.isinside(p):
                        fl.remove(f)
                        found = True
                        ans.append(f)
                        break
                if not found:
                    ans.append(None)

        else:
            for p in pl:
                found = False
                for f in self.faces:
                    if f.external:
                        continue
                    if f.isinside(p):
                        found = True
                        ans.append(f)
                        break
                if not found:
                    ans.append(None)

        return ans

    def load(self, filename):
        """reads a dcel from file using xygraph file format"""
        a = Xygraph.load(self, filename)
        self.build_dcel()
        return a

    def areas(self):
        return [f.area() for f in self.faces]

    def perimeters(self):
        return [f.perimeter() for f in self.faces]

    def nfaces(self):
        return len(self.faces)

    def nvertices(self):
        return len(self.vertices)

    def nedges(self):
        return len(self.hedges) / 2

    def delete_edge(self, v1, v2):
        """
        delete an edge (2 twin half-edges) from the dcel and update the whole dcel to be sane.
        :param v1: the first vertex of an edge to be deleted
        :param v2: the second vertex of an edge to be deleted
        :return: NA
        """
        # find a half-edge that joins these two points by iterating over the dcel
        rel_hedge = None
        for v1_hedge in v1.hedgelist:
            if v1_hedge.twin.origin == v2:
                rel_hedge = v1_hedge
        # find the faces adjacent to these half-edges
        face_1 = rel_hedge.face
        face_2 = rel_hedge.twin.face
        # remove connections to and from the half-edges, and connect neighbours to each other.
        rel_hedge.nexthedge.prevhedge = rel_hedge.twin.prevhedge
        rel_hedge.twin.prevhedge.nexthedge = rel_hedge.nexthedge
        rel_hedge.prevhedge.nexthedge = rel_hedge.twin.nexthedge
        rel_hedge.twin.nexthedge.prevhedge = rel_hedge.prevhedge
        # find a hedge on the boundary of the new face that was formed from merging prev faces.
        boundary_hedge = rel_hedge.prevhedge
        # delete half-edges from dcel structure and hedgelist of both vertices
        v1_vl_index = self.vl.index([v1.x, v1.y])
        v2_vl_index = self.vl.index([v2.x, v2.y])
        e_el_index = self.el.index([v1_vl_index, v2_vl_index]) if v1_vl_index < v2_vl_index else self.el.index(
            [v2_vl_index, v1_vl_index])
        v1.hedgelist.remove(rel_hedge)
        v2.hedgelist.remove(rel_hedge.twin)
        # sort the half-edges at v1 and v2 once again after deletion of half-edges
        v1.sortincident()
        v2.sortincident()
        self.hedges.remove(rel_hedge)
        self.hedges.remove(rel_hedge.twin)
        self.el.pop(e_el_index)
        # initialize a new face and put the half-edges in its boundary,i.e., change face pointers of half-edges.
        new_face = Face()
        new_face.wedge = boundary_hedge
        iter_hedge = boundary_hedge
        iter_hedge.face = new_face
        while iter_hedge.nexthedge != boundary_hedge:
            iter_hedge = iter_hedge.nexthedge
            iter_hedge.face = new_face
        # remove old faces and add new face
        self.faces.remove(face_1)
        self.faces.remove(face_2)
        self.faces.append(new_face)
        # complete everything related to that face.
        new_face.external = new_face.area() < 0

    def add_edge(self, v1, v2, legal=True):
        """
        add edge (2 half-edges) between two vertices v1 and v2 and change the dcel structure accordingly.
        :param v1: first vertex of the edge to be added
        :param v2: second vertex of the edge to be added
        :param legal: if the added edges are legal or not. defaults to true.
        :return: NA
        """
        # instantiate new half-edges
        h1 = Hedge(v1, v2)
        h2 = Hedge(v2, v1)
        h1.legal = legal
        h2.legal = legal
        h1.twin = h2
        h2.twin = h1

        # find an out half-edge from v1 that is on the boundary of the face that will get divided
        v1_out = None
        for v1_hedge in v1.hedgelist:
            iter_hedge = v1_hedge
            while iter_hedge.nexthedge.origin != v1:
                iter_hedge = iter_hedge.nexthedge
                if iter_hedge.external:
                    break
                if iter_hedge.origin == v2:
                    v1_out = v1_hedge
                    break
            if v1_out is not None:
                break
        # add half-edges to the half-edgelist of v1 and v2
        v1_in = v1_out.prevhedge
        v1.hedgelist.append(h2)
        v1.sortincident()
        v2.hedgelist.append(h1)
        v2.sortincident()

        # add half-edges to the list of half-edges
        self.hedges.append(h2)
        self.hedges.append(h1)

        # find an in half-edge incident on v2 that is on the boundary of the face that will get divided.
        v2_in = v1_out
        while v2_in.nexthedge.origin != v2:
            v2_in = v2_in.nexthedge

        # change the connections and add connections to new half-edges
        v2_out = v2_in.nexthedge

        v2_in.nexthedge = h1
        h1.prevhedge = v2_in
        h1.nexthedge = v1_out
        v1_out.prevhedge = h1

        v2_out.prevhedge = h2
        h2.nexthedge = v2_out
        v1_in.nexthedge = h2
        h2.prevhedge = v1_in
        # remove the original face from faces list
        self.faces.remove(v1_out.face)
        # instantiate new faces and put half-edges on its boundary according to the previous pointers to half-edges.
        new_face1 = Face()
        new_face1.wedge = v1_out
        iter_hedge = v1_out
        iter_hedge.face = new_face1
        while iter_hedge.nexthedge != v1_out:
            iter_hedge = iter_hedge.nexthedge
            iter_hedge.face = new_face1

        new_face2 = Face()
        new_face2.wedge = v1_in
        iter_hedge = v1_in
        iter_hedge.face = new_face2
        while iter_hedge.nexthedge != v1_in:
            iter_hedge = iter_hedge.nexthedge
            iter_hedge.face = new_face2

        new_face1.external = new_face1.area() < 0
        new_face2.external = new_face2.area() < 0
        # append new faces to the face list.
        self.faces.append(new_face1)
        self.faces.append(new_face2)
        # change the edge list for plotting
        v1_vl_index = self.vl.index([v1.x, v1.y])
        v2_vl_index = self.vl.index([v2.x, v2.y])
        self.el.append([v1_vl_index, v2_vl_index])
        # check dcel for dangling hedges.
        checkhedges(self.hedges)

    def incircle(self, v1, v2, v3, q):
        """
        tests if q is in/on/outside the cirle formed by v1, v2, v3 in that order.
        :param v1: a vertex object
        :param v2: a vertex object
        :param v3: a vertex object
        :param q: a query vertex object.
        v1, v2, v3 must be in counter-clockwise order.
        :return:
            1 if q is inside the circle
            0 if q is on the circle
            -1 if q is outside the circle
        """
        # v1, v2, v3 in that order are counter-clockwise in direction.
        incircle_det = np.array([[v1.x, v1.y, (v1.x ** 2) + (v1.y ** 2), 1],
                                 [v2.x, v2.y, (v2.x ** 2) + (v2.y ** 2), 1],
                                 [v3.x, v3.y, (v3.x ** 2) + (v3.y ** 2), 1],
                                 [q.x, q.y, (q.x ** 2) + (q.y ** 2), 1]], dtype=np.float32)

        det_val = np.linalg.det(incircle_det)
        if det_val > 0:
            return 1
        elif det_val == 0:
            return 0
        else:
            return -1

    def diagonal_triangulation_check(self, vj, vk, popped):
        """
        check if a diagonal can be added during the monotone polygon triangulation algorithm
        :param vj: vertex that is being checked
        :param vk: vertex that is on the stack
        :param popped: the vertex that was just popped from the stack
        :return:
            True if a diagonal can be added
            False if a diagonal can't be added

        Essentially checking if the edge connecting the popped vertex to vj intersects the edge connecting vj and vk.
        """
        ret = (lefton(popped, vj, vk) == lefton(vk, popped, vj)) and self.diagonal_fie(vk, vj)
        return ret

    @staticmethod
    def diagonal_fie(v1, v2):
        """
        The in-cone test
        :param v1: first vertex on the polygon
        :param v2: second vertex on the polygon
        :return:
            True if the diagonal between them lies incone
            False if the diagonal between them lies out of cone.
        """
        v1_ccw = None
        for v1_hedge in v1.hedgelist:
            iter_hedge = v1_hedge
            while iter_hedge.nexthedge.origin != v1:
                iter_hedge = iter_hedge.nexthedge
                if iter_hedge.external:
                    break
                if iter_hedge.origin == v2:
                    v1_ccw = v1_hedge
                    break
            if v1_ccw is not None:
                break

        v1_prev = v1_ccw.prevhedge
        if lefton(v1, v1_ccw.twin.origin, v1_prev.origin):
            return left(v1, v2, v1_prev.origin) and left(v2, v1, v1_ccw.twin.origin)
        else:
            return not (lefton(v1, v2, v1_ccw.twin.origin) and lefton(v2, v1, v1_prev.origin))

    def monotone_polygon_triangulation(self):
        """
        Does the monotone polygon triangulation assuming that the polygon is x-monotone.
        Assuming that there exists just one polygonal face in the dcel.
        :return:
        """
        # finding the interior face
        interior_face = None
        for i in self.faces:
            if not i.external:
                interior_face = i
        # finding the left and right extreme points
        max_x_vert = self.vertices[0]
        min_x_vert = self.vertices[0]
        for vert in self.vertices:
            if vert.x >= max_x_vert.x:
                if vert.x == max_x_vert.x and vert.y < max_x_vert.y:
                    max_x_vert = vert
                elif vert.x > max_x_vert.x:
                    max_x_vert = vert
            if vert.x <= min_x_vert.x:
                if vert.x == min_x_vert.x and vert.y < min_x_vert.y:
                    min_x_vert = vert
                elif vert.x < min_x_vert.x:
                    min_x_vert = vert
        # finding a half-edge that lies on the face boundary
        iter_hedge = None
        for hedge in max_x_vert.hedgelist:
            if hedge.face == interior_face:
                iter_hedge = hedge
        # finding the upper chain and lower chain of the face
        upper_chain = []
        lower_chain = []
        lower_chain_start = False
        upper_chain_start = False
        while iter_hedge.nexthedge.origin != max_x_vert:
            if iter_hedge.origin == max_x_vert:
                lower_chain_start = True
                upper_chain_start = False
                if len(upper_chain) > 0:
                    upper_chain.append(iter_hedge.origin)
                iter_hedge = iter_hedge.nexthedge
                continue
            if iter_hedge.origin == min_x_vert:
                upper_chain_start = True
                upper_chain.append(iter_hedge.origin)
                lower_chain_start = False
                iter_hedge = iter_hedge.nexthedge
                continue
            if lower_chain_start:
                lower_chain.append(iter_hedge.origin)
                iter_hedge = iter_hedge.nexthedge
                continue
            if upper_chain_start:
                upper_chain.append(iter_hedge.origin)
                iter_hedge = iter_hedge.nexthedge
                continue
            iter_hedge = iter_hedge.nexthedge
        upper_chain.append(iter_hedge.origin)
        upper_chain.append(max_x_vert)
        lower_chain = lower_chain[::-1]

        # merging the sorted order of upper and lower chains for the actual algorithm
        upper_pointer = 0
        lower_pointer = 0
        merged_chain = []
        status_list = []  # true if from upper chain, false if from lower chain
        while upper_pointer < len(upper_chain) and lower_pointer < len(lower_chain):
            if upper_chain[upper_pointer].x <= lower_chain[lower_pointer].x:
                merged_chain.append(upper_chain[upper_pointer])
                status_list.append(True)
                upper_pointer += 1
                continue
            if upper_chain[upper_pointer].x > lower_chain[lower_pointer].x:
                merged_chain.append(lower_chain[lower_pointer])
                status_list.append(False)
                lower_pointer += 1
                continue
        while upper_pointer < len(upper_chain):
            merged_chain.append(upper_chain[upper_pointer])
            status_list.append(True)
            upper_pointer += 1
        while lower_pointer < len(lower_chain):
            merged_chain.append(lower_chain[lower_pointer])
            status_list.append(False)
            lower_pointer += 1

        merged_chain = merged_chain[::-1]
        status_list = status_list[::-1]

        # start the triangulation
        stack = []  # the vertices are pushed in a tuple along with their index in the merged chain.
        stack.append((0, merged_chain[0]))
        stack.append((1, merged_chain[1]))
        for j in range(2, len(merged_chain) - 1):
            if status_list[j] ^ status_list[stack[-1][0]]:  # if the vertices are on opposite chains
                while len(stack) != 1:
                    current_top = stack.pop(-1)
                    self.add_edge(merged_chain[j], current_top[1])
                stack.pop(-1)
                stack.append((j - 1, merged_chain[j - 1]))
                stack.append((j, merged_chain[j]))
            else:  # if the vertices are on the same chain
                current_top = stack.pop(-1)
                # keep checking until you find the vertex that can't form a valid diagonal with vj
                while self.diagonal_triangulation_check(merged_chain[j], stack[-1][1], current_top[1]):
                    current_top = stack.pop(-1)
                    self.add_edge(current_top[1], merged_chain[j])
                    if len(stack) == 0:
                        break
                stack.append(current_top)
                stack.append((j, merged_chain[j]))
        stack.pop(-1)
        # add edges form the left most point to all points left in the stack.
        while len(stack) != 1:
            current_top = stack.pop(-1)
            self.add_edge(merged_chain[-1], current_top[1])
        print("triangulation done")

    def lawson_delaunay(self):
        # import pdb;pdb.set_trace()
        """
        Does lawson delaunay triangulation of the randomly triangulated polygon
        :return:
        """
        num_hedges = len(self.hedges)
        iteration = 0
        iterator = -1
        # iterate until convergence
        while True:
            iteration += 1
            iterator += 1
            iterator %= num_hedges
            # iterate over half-edges
            hedge = self.hedges[iterator]
            # find a non-legal edge
            if hedge.external or hedge.legal:
                continue
            else:
                # implement the incircle test
                v3 = hedge.origin
                v1 = hedge.nexthedge.origin
                v2 = hedge.prevhedge.origin
                q = hedge.twin.prevhedge.origin
                # test legality
                if self.incircle(v1, v2, v3, q) == 1:
                    self.delete_edge(v1, v3)
                    self.add_edge(v2, q)
                else:
                    hedge.legal = True
            # self.saveplot(f"ex1_{iteration}.eps")
            # break testing
            count_hedges = 0
            for hedge in self.hedges:
                if hedge.legal == True:
                    count_hedges += 1
            if count_hedges == num_hedges:
                break
        checkhedges(self.hedges)
        print(f"Edge Swap Algorithm Finished in {iteration} iterations over")


# Misc. functions


def hsort(h1, h2):
    """Sorts two half edges counterclockwise"""

    if h1.angle < h2.angle:
        return -1
    elif h1.angle > h2.angle:
        return 1
    else:
        return 0


def checkhedges(hl):
    """Consistency check of a hedge list: nexthedge, prevhedge"""

    for h in hl:
        if h.nexthedge not in hl or h.prevhedge not in hl:
            raise DcelError("Problems with an orphan hedge...")


def area2(a0, a1, query):
    """Determines the area of the triangle formed by a hedge and
    an external point"""

    # pa = hedge.twin.origin
    # pb = hedge.origin
    # pc = point
    return (a1.x - a0.x) * (query.y - a0.y) - (query.x - a0.x) * (a1.y - a0.y)


def lefton(a0, a1, query):
    """Determines if a point is to the left of a hedge"""

    return area2(a0, a1, query) >= 0


def left(a0, a1, query):
    """Determines if a point is to the left of a hedge"""

    return area2(a0, a1, query) > 0


def hangle(dx, dy):
    """Determines the angle with respect to the x axis of a segment
    of coordinates dx and dy
    """

    l = m.sqrt(dx * dx + dy * dy)
    if dy > 0:
        return m.acos(dx / (l + 1e-7))
    else:
        return 2 * m.pi - m.acos(dx / (l + 1e-7))


if __name__ == '__main__':
    import sys

    d = Dcel()
    d.load(sys.argv[1])
    for a, p in zip(d.areas(), d.perimeters()):
        print(a, p)
