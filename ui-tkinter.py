from tkinter import *
import tkinter as tk
from tkinter import ttk
from scipy.spatial import ConvexHull
import numpy as np
from dcel.dcel import Dcel
import itertools
import copy


class GUI:
    """
    The GUI class
    """
    def __init__(self, win_size=700):
        """
        Initialize the class attributes
        :param win_size: size of the window for the gui

        :attr lower_hull_vertices: the vertices of the lower hull of the convex-hull of points. List of lists.
        :attr upper_hull_vertices: the vertices of the upper hull of the convex-hull of points. List of lists.
        :attr win: the TK gui window object
        :attr points: the list of points. List of lists.
        :attr edge_list: the list of edges between points.
        :attr canvas: the canvas object of the TK module that helps draw on ui.
        :attr dcel_face1: the upper x-monotone polygon DCEL.
        :attr dcel_face2: the lower x-monotone polygon DCEL.

        :attr dcel: the DCEL for the completely triangulated set of points. Used after merging
        dcel_face1 and dcel_face2.
        """
        self.lower_hull_vertices = None
        self.upper_hull_vertices = None
        self.win = Tk()
        self.win.geometry(f"{win_size}x{win_size}")
        self.points = []
        self.edge_list = []
        self.canvas = Canvas(self.win, width=win_size, height=win_size, background="white")
        self.canvas.grid(row=0, column=0)
        self.canvas.bind('<p>', self.partition)
        self.canvas.bind('<r>', self.random_triangulate)
        self.canvas.bind('<d>', self.delaunay_triangulate)
        self.canvas.focus_set()
        self.canvas.bind('<Button-1>', self.draw_dot)

        self.dcel_face1 = Dcel()
        self.dcel_face2 = Dcel()

        self.dcel = Dcel()

        self.win.mainloop()

    def draw_dot(self, event):
        """
        Helper function to draw points and label them.
        :param event: mouse button press event.
        :return: NA
        """
        x1 = event.x
        y1 = event.y
        self.points.append([x1, y1])
        # Draw an oval in the given co-ordinates
        self.canvas.delete('all')
        for point in self.points:
            self.canvas.create_oval(point[0], point[1], point[0], point[1], fill="black", width=5)

    def draw_lines(self, vl, el, colour="black"):
        """
        helper function to draw edges between points.
        :param vl: list of vertices
        :param el: list of edges. list of list, containing the indices of vertices from the above list
        that each edge connects.

        :param colour: the colour of the edge to be plotted.
        :return: NA
        """
        for e in el:
            self.canvas.create_line(vl[e[0]][0], vl[e[0]][1], vl[e[1]][0], vl[e[1]][1], fill=colour)

    def partition(self, event):
        """
        Find the convex hull and partition the convex-polygon into two x-monotone polygons.
        Initialize their respective DCEL structures.
        :param event: keyboard press event.
        :return: NA
        """
        # remove duplicate points
        self.points.sort()
        self.points = list(self.points for self.points, _ in itertools.groupby(self.points))
        # do delaunay triangulation
        for point in self.points:
            self.canvas.create_text(point[0], point[1], text=f"{point[0]}, {point[1]}")
            self.canvas.create_oval(point[0], point[1], point[0], point[1], fill="black", width=5)
        # find convex hull
        hull = ConvexHull(np.array(self.points))
        points_in_order = np.array(self.points)[hull.vertices]
        # get x-extreme points of convex hull boundary
        ch_right_index = np.argmax(points_in_order, axis=0)[0]
        right_most_index = self.points.index(points_in_order[ch_right_index].tolist())
        ch_left_index = np.argmin(points_in_order, axis=0)[0]
        left_most_index = self.points.index(points_in_order[ch_left_index].tolist())

        # find upper and lower hull vertices
        upper_hull_vertices = []
        lower_hull_vertices = []
        upper_hull_flag = False
        lower_hull_flag = False
        for i in list(range(len(points_in_order))) * 2:
            if i == ch_right_index:
                lower_hull_flag = True
                upper_hull_flag = False
                if points_in_order[i].tolist() not in lower_hull_vertices:
                    lower_hull_vertices.append(points_in_order[i].tolist())
                if len(upper_hull_vertices) > 0 and points_in_order[i].tolist() not in upper_hull_vertices:
                    upper_hull_vertices.append(points_in_order[i].tolist())
                continue
            if upper_hull_flag and len(upper_hull_vertices) + len(lower_hull_vertices) < len(points_in_order) + 1:
                upper_hull_vertices.append(points_in_order[i].tolist())
                continue
            if i == ch_left_index:
                lower_hull_flag = False
                upper_hull_flag = True
                if len(lower_hull_vertices) > 0 and points_in_order[i].tolist() not in lower_hull_vertices:
                    lower_hull_vertices.append(points_in_order[i].tolist())
                if points_in_order[i].tolist() not in upper_hull_vertices:
                    upper_hull_vertices.append(points_in_order[i].tolist())
                continue
            if lower_hull_flag and len(upper_hull_vertices) + len(lower_hull_vertices) < len(points_in_order) + 1:
                lower_hull_vertices.append(points_in_order[i].tolist())
                continue

        # find the set of internal points
        internal_points = copy.deepcopy(self.points)
        for i in range(len(points_in_order)):
            internal_points.remove(points_in_order[i].tolist())
        # sort and find the order of points by x-coordinate first and y-coordinate second
        internal_points.sort(key=lambda x: (x[0], x[1]))

        # add internal points to upper and lower hull
        for internal_point in internal_points[::-1]:
            upper_hull_vertices.append(internal_point)
        for internal_point in internal_points:
            lower_hull_vertices.append(internal_point)

        self.upper_hull_vertices = upper_hull_vertices
        self.lower_hull_vertices = lower_hull_vertices
        self.dcel_face1.vl = self.upper_hull_vertices
        self.dcel_face1.el = [[i, i + 1] for i in range(len(self.upper_hull_vertices) - 1)] + [
            [len(self.upper_hull_vertices) - 1, 0]]
        self.dcel_face1.build_dcel()
        print("upper face built")
        self.canvas.delete('all')

        # plotting
        for point in self.points:
            self.canvas.create_text(point[0], point[1], text=f"{point[0]}, {point[1]}")
            self.canvas.create_oval(point[0], point[1], point[0], point[1], fill="black", width=5)
        self.draw_lines(self.dcel_face1.vl, self.dcel_face1.el, "red")
        self.dcel_face2.vl = self.lower_hull_vertices
        self.dcel_face2.el = [[i, i + 1] for i in range(len(self.lower_hull_vertices) - 1)] + [
            [len(self.lower_hull_vertices) - 1, 0]]
        self.dcel_face2.build_dcel()
        self.draw_lines(self.dcel_face2.vl, self.dcel_face2.el, "blue")
        print("lower face built")

        # ui essentials
        self.canvas.unbind('<Button-1>')

    def random_triangulate(self, event):
        """
        randomly triangulate both the polygons using monotone polygon triangulation. Then merge the 2 DCEL data
         structures and initialize the one DCEL data structure for all points
        :param event: keyboard press event.
        :return: NA
        """
        self.dcel_face1.monotone_polygon_triangulation()
        self.dcel_face2.monotone_polygon_triangulation()

        first_vl = self.dcel_face1.vl
        second_vl = self.dcel_face2.vl
        first_el = self.dcel_face1.el
        second_el = self.dcel_face2.el

        final_vl = copy.deepcopy(first_vl)
        for v in second_vl:
            if v not in final_vl:
                final_vl.append(v)
        final_el = copy.deepcopy(first_el)
        for e in second_el:
            new_e = [final_vl.index(second_vl[e[0]]), final_vl.index(second_vl[e[1]])]
            final_el.append(new_e)
        final_el = [e[::-1] if e[0]>e[1] else e for e in final_el]
        final_el.sort(key=lambda x:(x[0], x[1]))
        res = []
        [res.append(x) for x in final_el if x not in res]
        final_el = res
        self.dcel.vl = final_vl
        self.dcel.el = final_el
        self.dcel.minmax()
        self.dcel.build_dcel()
        print("triangulations merged")

        # plotting
        self.canvas.delete('all')
        for point in self.points:
            self.canvas.create_text(point[0], point[1], text=f"{point[0]}, {point[1]}")
            self.canvas.create_oval(point[0], point[1], point[0], point[1], fill="black", width=5)
        self.draw_lines(self.dcel.vl, self.dcel.el, colour="blue")

        # ui essentials
        self.canvas.unbind('<p>')

    def delaunay_triangulate(self, event):
        """
        Perform lawson swap on the complete DCEL structure that is randomly triangulated.
        :param event: keyboard press event
        :return: NA
        """
        self.dcel.lawson_delaunay()
        # plotting
        self.canvas.delete('all')
        for point in self.points:
            self.canvas.create_text(point[0], point[1], text=f"{point[0]}, {point[1]}")
            self.canvas.create_oval(point[0], point[1], point[0], point[1], fill="black", width=5)
        self.draw_lines(self.dcel.vl, self.dcel.el, colour="green")

        # ui essentials
        self.canvas.unbind('<r>')

# Create a canvas widget
gui = GUI()
