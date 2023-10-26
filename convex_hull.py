from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT6':
    from PyQt6.QtCore import QLineF, QPointF, QObject
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time

# Some global color constants that might be useful
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)

# Global variable that controls the speed of the recursion automation, in seconds
#
PAUSE = 0.25


# This function takes in the hulls, finds the upper and lower tangent and re-points the new tangent endpoints
# to the opposite hull and throws out the middle points that are no longer necessary. n timei   in
def merge(left_half, right_half):
    # Prepare lower and upper tangents
    upper_left = max(left_half, key=lambda p: p.point.x())
    upper_right = min(right_half, key=lambda p: p.point.x())
    bottom_left = upper_left
    bottom_right = upper_right

    # Finds upper tangent and stores points in upper_left and upper_right
    while True:
        last_l = upper_left
        last_r = upper_right
        if upper_right.rotate_cw:
            while direction(upper_left, upper_right, upper_right.rotate_cw) < 0:
                upper_right = upper_right.rotate_cw
        if upper_left.rotate_ccw:
            while direction(upper_right, upper_left, upper_left.rotate_ccw) > 0:
                upper_left = upper_left.rotate_ccw

        if upper_left == last_l and upper_right == last_r:
            break

    # Finds lower tangent and stores points in lower_left and lower_right
    while True:
        last_l = bottom_left
        last_r = bottom_right
        if bottom_right.rotate_ccw:
            while direction(bottom_left, bottom_right, bottom_right.rotate_ccw) > 0:
                bottom_right = bottom_right.rotate_ccw
        if bottom_left.rotate_cw:
            while direction(bottom_right, bottom_left, bottom_left.rotate_cw) < 0:
                bottom_left = bottom_left.rotate_cw
        if bottom_right == last_r and bottom_left == last_l:
            break

    # Adjusts the pointers to the correct tangent points to eliminate the middle points
    upper_left.rotate_cw = upper_right
    upper_right.rotate_ccw = upper_left
    bottom_left.rotate_ccw = bottom_right
    bottom_right.rotate_cw = bottom_left

    # Creates array of the new combined hull's points
    start = upper_left
    result = []
    while True:
        result.append(upper_left)
        upper_left = upper_left.rotate_ccw
        if upper_left == start:
            break
    return result


class ConvexHullSolver(QObject):

    # Class constructor
    def __init__(self):
        super().__init__()
        self.pause = False

    # Some helper methods that make calls to the GUI, allowing us to send updates
    # to be displayed.

    def showTangent(self, line, color):
        self.view.addLines(line, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseTangent(self, line):
        self.view.clearLines(line)

    def blinkTangent(self, line, color):
        self.showTangent(line, color)
        self.eraseTangent(line)

    def showHull(self, polygon, color):
        self.view.addLines(polygon, color)
        if self.pause:
            time.sleep(PAUSE)

    def eraseHull(self, polygon):
        self.view.clearLines(polygon)

    def showText(self, text):
        self.view.displayStatusText(text)

    # This is the method that gets called by the GUI and actually executes
    # the finding of the hull
    def compute_hull(self, points, pause, view):
        self.pause = pause
        self.view = view
        assert (type(points) == list and type(points[0]) == QPointF)

        t1 = time.time()
        # TODO: SORT THE POINTS BY INCREASING X-VALUE
        points = sorted(points, key=lambda p: p.x())
        t2 = time.time()

        t3 = time.time()
        # this is a dummy polygon of the first 3 unsorted points
        # polygon = [QLineF(points[i], points[(i + 1) % 3]) for i in range(3)]
        convex_hull = self.convex_hull(points)
        # TODO: REPLACE THE LINE ABOVE WITH A CALL TO YOUR DIVIDE-AND-CONQUER CONVEX HULL SOLVER
        convex_hull_points = []
        for point in convex_hull:
            convex_hull_points.append(point.point)
        polygon = [QLineF(convex_hull_points[i], convex_hull_points[(i + 1)]) for i in
                   range(len(convex_hull_points) - 1)]
        polygon.append(QLineF(convex_hull_points[-1], convex_hull_points[0]))
        t4 = time.time()

        # when passing lines to the display, pass a list of QLineF objects.  Each QLineF
        # object can be created with two QPointF objects corresponding to the endpoints
        self.showHull(polygon, RED)
        self.showText('Time Elapsed (Convex Hull): {:3.3f} sec'.format(t4 - t3))

    # Splits an array of points down to 1 point and initiates the merge function. logn time
    def convex_hull(self, points):
        if len(points) == 1:
            points[0] = Point(points[0])
            return points
        # TODO: RIGHT HERE
        median = int(len(points) / 2)
        left_half = self.convex_hull(points[0: median])
        right_half = self.convex_hull(points[median:])
        return merge(left_half, right_half)


# Point class that store a QPointF object and pointers to the points on either side of the point. Constant Time
class Point:
    def __init__(self, p):
        self.point = p
        self.rotate_cw = None
        self.rotate_ccw = None

    def subtract(self, p):
        check = QPointF(self.point.x() - p.point.x(), self.point.y() - p.point.y())
        return check


# Creates vector to show slope direction. Constant Time
def cross_product(p1, p2):
    return p1.x() * p2.y() - p2.x() * p1.y()


# Finds a slope increase or decrease. Constant Time
def direction(p1, p2, p3):
    return cross_product(p3.subtract(p1), p2.subtract(p1))

