from matplotlib.artist import Artist
from numpy import sqrt, nonzero, equal, array, asarray, amin, searchsorted
import numpy as np

from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

def dist_point_to_segment(p, s0, s1):
    p = np.asarray(p, dtype=float)
    s0 = np.asarray(s0, dtype=float)
    s1 = np.asarray(s1, dtype=float)
    v = s1 - s0
    w = p - s0
    c1 = np.dot(w, v)
    if c1 <= 0:
        return np.hypot(p[0] - s0[0], p[1] - s0[1])
    c2 = np.dot(v, v)
    if c2 <= c1:
        return np.hypot(p[0] - s1[0], p[1] - s1[1])
    b = c1 / c2
    pb = s0 + b * v
    return np.hypot(p[0] - pb[0], p[1] - pb[1])

class PolygonInteractor(object):
    showverts = True
    epsilon = 5  # max distance to count as a vertex hit

    def __init__(self, ax_edit, poly, no_nodes, digitized_file, reso, Z, topo, out):
        if poly.figure is None:
            raise RuntimeError('You must first add the polygon to a figure or canvas before defining the interactor')
        self.ax_edit = ax_edit
        self.out_file = out
        canvas = poly.figure.canvas
        self.topo = topo
        self.poly = poly
        self.l = []
        self.b = []
        self.reso = reso
        self.Z = Z
        self.cum_l = []
        self.no_nodes = no_nodes
        self.cum_l = np.ndarray.tolist(np.cumsum(self.no_nodes))
        self.cum_l.insert(0, 0)
        
        # ===================================================================
        # Display separate closed boundaries without cross-layer connecting lines.
        # ===================================================================
        # Hide the patch outline to avoid connecting unrelated contours.
        self.poly.set_edgecolor('none')
        
        xy = self.poly.get_xy()
        x, y = xy[:, 0], xy[:, 1]
        self.locked_vertices = self._top_surface_vertices(xy)
        
        # Show draggable vertices without a connecting line.
        self.line = Line2D(x, y, marker='o', markersize=4, markerfacecolor='r', linestyle='none', zorder=10)
        self.ax_edit.add_line(self.line)
        
        # Add a separate display line behind editable control points.
        line_color = 'blue' if 'normal' in self.out_file else 'red'
        self.display_line = Line2D([], [], color=line_color, linewidth=1.5, zorder=1)
        self.ax_edit.add_line(self.display_line)
        
        # Initialize the closed, segmented display path.
        self.update_display_line()
        # ===================================================================

        self.cid = self.poly.add_callback(self.poly_changed)
        self._ind = None 

        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('button_press_event', self.save_poly)
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('button_release_event', self.button_release_callback)
        canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.canvas = canvas

    def _top_surface_vertices(self, xy):
        try:
            topo = np.loadtxt(self.topo)
            if topo.ndim == 1:
                topo = topo.reshape(1, -1)
            topo_x = topo[:, 0]
            topo_y = topo[:, 1] / 1000.0
        except Exception:
            return {}

        locked = {}
        x_tol = max(abs(float(self.reso)) * 0.51, 1.0e-9)
        y_tol = 0.05
        for i, (x, y) in enumerate(xy):
            nearest = int(np.argmin(np.abs(topo_x - x)))
            if abs(topo_x[nearest] - x) <= x_tol and abs(topo_y[nearest] - y) <= y_tol:
                locked[i] = (float(topo_x[nearest]), float(topo_y[nearest]))
        return locked

    def _boundary_x(self, x, original_x=None):
        xmin, xmax = sorted(self.ax_edit.get_xlim())
        tol = max(abs(float(self.reso)) * 0.5, 1.0e-9)

        if original_x is not None:
            if abs(original_x - xmin) <= tol:
                return xmin
            if abs(original_x - xmax) <= tol:
                return xmax

        if abs(x - xmin) <= tol:
            return xmin
        if abs(x - xmax) <= tol:
            return xmax
        return x

    def _snap_boundary_vertices(self, xy):
        snapped = np.array(xy, copy=True)
        for i in range(len(snapped)):
            snapped[i, 0] = self._boundary_x(snapped[i, 0])
            if i in self.locked_vertices:
                snapped[i] = self.locked_vertices[i]
        return snapped

    def update_display_line(self):
        # Separate layer contours with NaN points and close each contour.
        xy = self.poly.get_xy()
        x, y = xy[:, 0], xy[:, 1]
        
        disp_x = []
        disp_y = []
        
        for i in range(len(self.no_nodes)):
            start = self.cum_l[i]
            end = self.cum_l[i+1]
            
            layer_x = x[start:end].tolist()
            layer_y = y[start:end].tolist()
            
            disp_x.extend(layer_x)
            disp_y.extend(layer_y)
            
            if len(layer_x) > 0:
                # Close the current contour by repeating its first vertex.
                disp_x.append(layer_x[0])
                disp_y.append(layer_y[0])
                
            # Insert a NaN break so separate contours are not connected.
            disp_x.append(np.nan)
            disp_y.append(np.nan)
            
        self.display_line.set_data(disp_x, disp_y)

    def draw_callback(self, event):
        pass

    def poly_changed(self, poly):
        vis = self.line.get_visible()
        Artist.update_from(self.line, self.poly)
        self.line.set_visible(vis)  
        self.line.set_linestyle('none') # Disable the control-point line.
        self.update_display_line()      # Refresh the closed display contours.

    def save_poly(self, event):
        if event.button == 3:
            xy = self._snap_boundary_vertices(self.poly.get_xy())
            self.poly.set_xy(xy)
            self.line.set_data(xy[:, 0], xy[:, 1])
            self.update_display_line()
            self.l = xy[:, 0]
            self.b = xy[:, 1]
            print("total length:", len(self.l))
            
            with open(self.out_file, "w") as f:
                f.write(">  \t  %d\n " % (len(self.no_nodes)))
                for i in range(len(self.no_nodes)):
                    f.write(">  \t  %d\n " % ((self.no_nodes[i])))
                    for j in range(self.cum_l[i], self.cum_l[i+1]):
                        f.write("%f \t %f \n" % (self.l[j], self.b[j]))

    def get_ind_under_point(self, event):
        xy = np.asarray(self.poly.get_xy())
        xyt = self.poly.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = np.hypot(xt - event.x, yt - event.y)
        
        min_d = np.amin(d)
        
        # Identify coincident or nearly coincident vertices.
        indseq = np.nonzero(d - min_d < 0.1)[0]
        indseq = np.array([idx for idx in indseq if idx not in self.locked_vertices], dtype=int)
        
        ind = indseq
        print("Grabbed overlapping point indices:", ind)
        if len(ind) == 0:
            return None
        if d[ind[0]] >= self.epsilon:
            ind = None
        return ind

    def button_press_callback(self, event):
        if not self.showverts: return
        if event.inaxes is None: return
        if event.button != 1: return
        self._ind = self.get_ind_under_point(event)

    def button_release_callback(self, event):
        if not self.showverts: return
        if event.button != 1: return
        self._ind = None

    def key_press_callback(self, event):
        if not event.inaxes: return
        if event.key == 't':
            self.showverts = not self.showverts
            self.line.set_visible(self.showverts)
            if not self.showverts: self._ind = None
            
        elif event.key == 'd':
            ind = self.get_ind_under_point(event)
            if ind is not None:
                sorted_indices = sorted(ind, reverse=True)
                xy_list = list(self.poly.get_xy())
                
                for idx in sorted_indices:
                    if idx < len(xy_list):
                        xy_list.pop(idx)
                        here = searchsorted(self.cum_l, idx)
                        for jj in range(here, len(self.cum_l)):
                            self.cum_l[jj] -= 1
                            
                self.poly.set_xy(xy_list)
                new_xy = self.poly.get_xy()
                self.line.set_data(new_xy[:, 0], new_xy[:, 1])
                self.update_display_line()
                self.canvas.draw_idle()
                
        elif event.key == 'i':
            ind = self.get_ind_under_point(event)
            if ind is not None:
                sorted_indices = sorted(ind, reverse=True)
                xy_list = list(self.poly.get_xy())
                
                for idx in sorted_indices:
                    if idx > 0:
                        curr_x, curr_y = xy_list[idx]
                        new_x = curr_x - self.reso
                        new_y = curr_y
                        
                        xy_list.insert(idx, (new_x, new_y))
                        
                        here = searchsorted(self.cum_l, idx)
                        for jj in range(here, len(self.cum_l)):
                            self.cum_l[jj] += 1
                            
                self.poly.set_xy(xy_list)
                new_xy = self.poly.get_xy()
                self.line.set_data(new_xy[:, 0], new_xy[:, 1])
                self.update_display_line()
                self.canvas.draw_idle()
                
        return

    def motion_notify_callback(self, event):
        if not self.showverts: return
        if self._ind is None: return
        if event.inaxes is None: return
        if event.button != 1: return
        
        # x = round(event.xdata / float(self.reso)) * float(self.reso)
        # y = event.ydata
        
        # xy = self.poly.get_xy()
        
        # Keep the original logic for moving coincident vertices together.
        # for idx in self._ind:
        #     xy[idx] = [x, y]
            
        # self.poly.set_xy(xy)
        # self.line.set_data(xy[:, 0], xy[:, 1])
        
        # Update the display line after moving vertices.
        # self.update_display_line()
        
        # self.canvas.draw_idle()

        # Use the raw mouse coordinates for editing.
        xy = self.poly.get_xy()
        original_x = xy[self._ind[0], 0] if len(self._ind) > 0 else None
        x = self._boundary_x(event.xdata, original_x)
        y = event.ydata
        
        for idx in self._ind:
            if idx in self.locked_vertices:
                continue
            xy[idx] = [x, y]
            
        self.poly.set_xy(xy)
        self.line.set_data(xy[:, 0], xy[:, 1])
        self.update_display_line()
        self.canvas.draw_idle()
