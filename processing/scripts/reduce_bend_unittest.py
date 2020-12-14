"""
Unit test for reduce_bend algorithm
"""

import unittest
from qgis.core import QgsApplication
from algo_reduce_bend import reduce_bends, is_polygon, is_line_string
from qgis.core import QgsPoint, QgsLineString, QgsPolygon, QgsFeature, QgsGeometry

def qgs_line_string_to_xy(qgs_line_string):

    qgs_points = qgs_line_string.points()
    lst_x = []
    lst_y = []
    for qgs_point in qgs_points:
        lst_x.append(qgs_point.x())
        lst_y.append(qgs_point.y())

    return (lst_x, lst_y)

def plot_lines(qgs_line_string, qgs_new_line):

    line0_lst_x, line0_lst_y = qgs_line_string_to_xy(qgs_line_string)
    line1_lst_x, line1_lst_y = qgs_line_string_to_xy(qgs_new_line)

    import matplotlib.pyplot as plt
    plt.plot(line0_lst_x, line0_lst_y, 'b')
    plt.plot(line1_lst_x, line1_lst_y, 'r')
    plt.show()


def build_and_launch(title, qgs_geoms, diameter_tol, del_pol=False, del_hole=False):

    print(title)
    qgs_features = []
    for qgs_geom in qgs_geoms:
        qgs_feature = QgsFeature()
        qgs_feature.setGeometry(qgs_geom)
        qgs_features.append(qgs_feature)

    rb_results = reduce_bends(qgs_features, diameter_tol, None, del_pol, del_hole)
    qgs_features_out = rb_results.qgs_features_out

    qgs_geoms_out = []
    for qgs_feature_out in qgs_features_out:
        qgs_geoms_out.append(qgs_feature_out.geometry())

    return qgs_geoms_out

def create_line(coords, ret_geom=True):

    qgs_points = []
    for coord in coords:
        qgs_points.append(create_point(coord, False))

    if ret_geom:
        ret_val = QgsGeometry(QgsLineString(qgs_points))
    else:
        ret_val = QgsLineString(qgs_points).clone()

    return ret_val

def create_point(coord, ret_geom=True):

    qgs_point = QgsPoint(coord[0], coord[1])
    if ret_geom:
        ret_val = QgsGeometry(qgs_point)
    else:
        ret_val = qgs_point.clone()

    return ret_val

def create_polygon(outer, inners):

#    outer_line = create_line(outer, False)
#    inner_lines = [create_line(inner, False) for inner in inners]
#    qgs_pol = QgsPolygon(outer_line, inner_lines)
#    qgs_geom = QgsGeometry(qgs_pol)

#    outer_line = QgsLineString([QgsPoint(0,0), QgsPoint(0,10), QgsPoint(10,10), QgsPoint(10,0), QgsPoint(0,0)])
#    inner_lines = [QgsLineString([QgsPoint(5,5), QgsPoint(5,6), QgsPoint(6,6), QgsPoint(6,5), QgsPoint(5,5)])]
#    qgs_geom = QgsGeometry(QgsPolygon(outer_line, inner_lines))




#    outer_line = QgsLineString([QgsPoint(0, 0), QgsPoint(0, 10), QgsPoint(10, 10), QgsPoint(10, 0), QgsPoint(0, 0)])
#    inner_lines = [QgsLineString([QgsPoint(5, 5), QgsPoint(5, 6), QgsPoint(6, 6), QgsPoint(6, 5), QgsPoint(5, 5)])]
#    qgs_pol = QgsPolygon(outer_line, inner_lines)
#    qgs_geom = QgsGeometry(qgs_pol.clone())

#    outer_line = QgsLineString([QgsPoint(0, 0), QgsPoint(0, 10), QgsPoint(10, 10), QgsPoint(10, 0), QgsPoint(0, 0)])
#    inner_lines = [QgsLineString([QgsPoint(5, 5), QgsPoint(5, 6), QgsPoint(6, 6), QgsPoint(6, 5), QgsPoint(5, 5)])]
    outer_line = create_line(outer, False)
    qgs_pol = QgsPolygon()
    qgs_pol.setExteriorRing(outer_line)
    for inner in inners:
        inner_line = create_line(inner, False)
        qgs_pol.addInteriorRing(inner_line)
    qgs_geom = QgsGeometry(qgs_pol)

    return qgs_geom


class Test(unittest.TestCase):
    """
    Class allowing to test the algorithm
    """

    def test_case01(self):
        title = "Test 01: Straight lines"
        qgs_geom0 = create_line([(0,0),(30,0)])
        qgs_geom1 = create_line([(0,10),(30,10)])
        qgs_geom2 = create_point((0,100))
        qgs_feature_out = build_and_launch(title,[qgs_geom0, qgs_geom1, qgs_geom2], 3)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        val1 = qgs_geom1.equals(qgs_feature_out[1])
        val2 = qgs_geom2.equals(qgs_feature_out[2])
        self.assertTrue (val0 and val1 and val2, title)

    def test_case02(self):

        title = "Test 02: Co-linear point"
        in_geom0 = create_line([(0, 0), (20, 0), (25.0000000001, 0.0000000001), (30, 0)])
        in_geom1 = create_line([(0, 10), (30, 10), (35.000000001, 10.00000000001), (40, 10)])
        in_geom2 = create_point((0, 100))
        out_geom0 = create_line([(0, 0), (30, 0)])
        out_geom1 = create_line([(0, 10), (40, 10)])
        qgs_feature_out = build_and_launch(title, [in_geom0, in_geom1, in_geom2], 3)
        val0 = out_geom0.equals(qgs_feature_out[0])
        val1 = out_geom1.equals(qgs_feature_out[1])
        val2 = in_geom2.equals(qgs_feature_out[2])
        self.assertTrue(val0 and val1 and val2, title)

    def test_case03(self):
        title = "Test 03: Small bend"
        in_geom0 = create_line([(0, 0), (30, 0)])
        in_geom1 = create_line([(0, 10), (30, 10), (30, 11), (31, 11), (31, 10), (40, 10), (50, 10), (50, 11), (51, 10), (60, 10)])
        in_geom2 = create_point((0, 100))
        out_geom0 = create_line([(0, 0), (30, 0)])
        out_geom1 = create_line([(0, 10), (60, 10)])
        qgs_feature_out = build_and_launch(title, [in_geom0, in_geom1, in_geom2], 3)
        val0 = out_geom0.equals(qgs_feature_out[0])
        val1 = out_geom1.equals(qgs_feature_out[1])
        val2 = in_geom2.equals(qgs_feature_out[2])
        self.assertTrue(val0 and val1 and val2, title)

    def test_case04(self):
        title = "Test 04: Polygon with bend"
        outer = [(0, 0), (0, 20), (10, 20), (10, 21), (11, 21), (11, 20), (20, 20), (20, 0), (10, -.1)]
        inner = [(5, 5), (5, 6), (6, 6), (6, 5)]
        in_geom0 = create_polygon(outer, [inner])
        outer = [(20,0), (10,-0.10000000000000001), (0,0), (0,20), (20,20), (20,0)]
        inner = [(5,5), (5,6), (6,6), (6,5), (5,5)]
        qgs_feature_out = build_and_launch(title, [in_geom0], 3)
        out_geom0 = create_polygon(outer, [inner])
        val0 = out_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case05(self):
        title = "Test 05: Polygon with line in bend"
        coord = [(0, 0), (0, 20), (10, 20), (10, 21), (11, 21), (11, 20), (20, 20), (20, 0)]
        qgs_geom0 = create_polygon(coord, [])
        qgs_geom1 = create_line([(10.1, 20.5), (10.2, 20.6), (10.3, 20.5)])
        qgs_feature_out = build_and_launch(title, [qgs_geom0, qgs_geom1], 3)
        coord = [(20,20), (20,0), (0,0), (0,20), (10,20), (10,21), (11,21), (11,20), (20,20)]
        out_geom0 = create_polygon(coord, [])
        out_geom1 = create_line([(10.1, 20.5), (10.3, 20.5)])
        val0 = out_geom0.equals(qgs_feature_out[0])
        val1 = out_geom1.equals(qgs_feature_out[1])
        self.assertTrue(val0 and val1, title)

    def test_case06(self):
        title = "Test 06: Polygon with point in bend"
        coord = [(0, 0), (0, 20), (10, 20), (10, 21), (11, 21), (11, 20), (20, 20), (20, 0)]
        qgs_geom0 = create_polygon(coord, [])
        qgs_geom1 = create_point((10.1,20.5))
        qgs_feature_out = build_and_launch(title, [qgs_geom0, qgs_geom1], 3)
        coord = [(20,20), (20,0), (0,0), (0,20), (10,20), (10,21), (11,21), (11,20), (20,20)]
        out_geom0 = create_polygon(coord, [])
        val0 = out_geom0.equals(qgs_feature_out[0])
        val1 = qgs_geom1.equals(qgs_feature_out[1])
        self.assertTrue(val0 and val1, title)

    def test_case07(self):
        title = "Test 07: Line String self intersecting after bend reduction"
        coord = [(0, 20), (10, 20), (10, 21), (11, 21), (11, 20), (30, 20), (30,0), (10.5,0), (10.5,20.5)]
        qgs_geom0 = create_line(coord)
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3)
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case08(self):
        title = "Test 08: Polygon with a bend containing a hole to delete"
        coord0 = [(0, 0), (0, 20), (10, 20), (10, 21), (11, 21), (11, 20), (20, 20), (20, 0)]
        coord1 = [(10.1, 20.1), (10.1, 20.2), (10.2, 20.2), (10.2, 20.1), (10.1, 20.1)]
        qgs_geom0 = create_polygon(coord0, [coord1])
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3, del_pol=True, del_hole=True)
        coord = [(20,20), (20,0), (0,0), (0,20), (20,20)]
        out_geom0 = create_polygon(coord, [])
        val0 = out_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case09(self):
        title = "Test 09: Small polygon no bend with one small hole. Hole is deleted"
        coord0 = [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]
        coord1 = [(0.1,0.1), (0.1,0.2), (0.2,0.2), (0.2,0.1), (0.1,0.1)]
        qgs_geom0 = create_polygon(coord0, [coord1])
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3, del_pol=False, del_hole=True)
        coord0 = [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]
        qgs_geom0 = create_polygon(coord0, [])
        val0 = qgs_geom0.equals(qgs_feature_out[0])
        self.assertTrue(val0, title)

    def test_case10(self):
        title = "Test 10: Small polygon no bend with one small hole. Feature is deleted"
        coord0 = [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]
        coord1 = [(0.1,0.1), (0.1,0.2), (0.2,0.2), (0.2,0.1), (0.1,0.1)]
        qgs_geom0 = create_polygon(coord0, [coord1])
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3, del_pol=True, del_hole=False)
        self.assertEqual(len(qgs_feature_out), 0, title)

    def test_case11(self):
        title = "Test 11: Small polygon no bend with one small hole. Feature is deleted"
        coord0 = [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]
        coord1 = [(0.1,0.1), (0.1,0.2), (0.2,0.2), (0.2,0.1), (0.1,0.1)]
        qgs_geom0 = create_polygon(coord0, [coord1])
        qgs_feature_out = build_and_launch(title, [qgs_geom0], 3, del_pol=True, del_hole=True)
        self.assertEqual(len(qgs_feature_out), 0, title)






# Supply path to qgis install location
QgsApplication.setPrefixPath("/usr/bin/qgis", False)

# profile_folder = 'C:\\Users\\berge\\AppData\\Roaming\\QGIS\\QGIS3\\profiles\\test12'
profile_folder = '.'
# Create a reference to the QgsApplication.  Setting the second argument to False disables the GUI.
app = QgsApplication([], False, profile_folder)

# Load providers
app.initQgis()