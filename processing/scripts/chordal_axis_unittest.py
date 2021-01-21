"""
Unit test for reduce_bend algorithm
"""

import unittest
from algo_chordal_axis import ChordalAxis, GenUtil
from qgis.core import QgsApplication
from qgis.core import QgsPoint, QgsLineString, QgsPolygon, QgsFeature, QgsGeometry, QgsProcessingFeedback

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
#    line1_lst_x, line1_lst_y = qgs_line_string_to_xy(qgs_new_line)

    import matplotlib.pyplot as plt
    plt.plot(line0_lst_x, line0_lst_y, 'b')
#    plt.plot(line1_lst_x, line1_lst_y, 'r')
    plt.show()


def build_and_launch(title, qgs_geoms, diameter_tol, del_pol=False, del_hole=False):

    print(title)
    qgs_features = []
    feedback = QgsProcessingFeedback()
    for qgs_geom in qgs_geoms:
        qgs_feature = QgsFeature()
        qgs_feature.setGeometry(qgs_geom)
        qgs_features.append(qgs_feature)

    rb_results = ReduceBend.reduce(qgs_features, diameter_tol, feedback, del_pol, del_hole, True)
    log = feedback.textLog()
    print (log)
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

    outer_line = create_line(outer, False)
    qgs_pol = QgsPolygon()
    qgs_pol.setExteriorRing(outer_line)
    for inner in inners:
        inner_line = create_line(inner, False)
        qgs_pol.addInteriorRing(inner_line)
    qgs_geom = QgsGeometry(qgs_pol)

    return qgs_geom

def coords_shift(pos, coords):

    new_coords = coords[pos:] + coords[0:pos-1] + [coords[pos]]

    return new_coords


class Test(unittest.TestCase):
    """
    Class allowing to test the algorithm
    """

    """
    def test_case00(self):
        f = open("/home/daneli/test.txt", "r")
        wkt = f.read()
        geom = QgsGeometry()
        geom1 = geom.fromWkt(wkt)
#        pol = geom.constGet().clone()
        title = "Test 01: 2 simple line segment, simple triangle  and one point"
#        qgs_geom0 = create_line([(0,0),(30,0)])
#        qgs_geom1 = create_polygon([(10,10),(15,20), (20,10), (10,10)], [])
#        qgs_geom2 = create_point((0,100))
        qgs_feature_out = build_and_launch(title,[geom1], 3)
#        val0 = qgs_geom0.equals(qgs_feature_out[0])
#        val1 = qgs_geom1.equals(qgs_feature_out[1])
#        val2 = qgs_geom2.equals(qgs_feature_out[2])
        self.assertTrue (True, title)

    """
    def test_case000_3(self):
        title = "Test 00: Empty file"
        txt_geom = "MultiPolygonZ (((559397.56532425712794065 5539480.39943350851535797 0, 559402.70997757744044065 5539469.33357657492160797 0, 559411.65791458915919065 5539480.49940909445285797 0, 559397.56532425712794065 5539480.39943350851535797 0)),((559411.65791458915919065 5539480.49940909445285797 0, 559402.70997757744044065 5539469.33357657492160797 0, 559417.69444413017481565 5539469.00062979757785797 0, 559411.65791458915919065 5539480.49940909445285797 0)),((559411.65791458915919065 5539480.49940909445285797 0, 559417.69444413017481565 5539469.00062979757785797 0, 559418.08039994072169065 5539469.02193106710910797 0, 559411.65791458915919065 5539480.49940909445285797 0)),((559426.15767044853419065 5539481.46608389914035797 0, 559411.65791458915919065 5539480.49940909445285797 0, 559418.08039994072169065 5539469.02193106710910797 0, 559426.15767044853419065 5539481.46608389914035797 0)),((559418.08039994072169065 5539469.02193106710910797 0, 559424.97224564384669065 5539469.94081534445285797 0, 559426.15767044853419065 5539481.46608389914035797 0, 559418.08039994072169065 5539469.02193106710910797 0)),((559426.15767044853419065 5539481.46608389914035797 0, 559424.97224564384669065 5539469.94081534445285797 0, 559429.37440628837794065 5539468.34004630148410797 0, 559426.15767044853419065 5539481.46608389914035797 0)),((559428.60399002861231565 5539482.51448477804660797 0, 559426.15767044853419065 5539481.46608389914035797 0, 559429.37440628837794065 5539468.34004630148410797 0, 559428.60399002861231565 5539482.51448477804660797 0)),((559428.60399002861231565 5539482.51448477804660797 0, 559429.37440628837794065 5539468.34004630148410797 0, 559443.25630326103419065 5539470.04579581320285797 0, 559428.60399002861231565 5539482.51448477804660797 0)),((559430.05711502861231565 5539484.61342276632785797 0, 559428.60399002861231565 5539482.51448477804660797 0, 559443.25630326103419065 5539470.04579581320285797 0, 559430.05711502861231565 5539484.61342276632785797 0)),((559441.87690872978419065 5539483.64760245382785797 0, 559430.05711502861231565 5539484.61342276632785797 0, 559443.25630326103419065 5539470.04579581320285797 0, 559441.87690872978419065 5539483.64760245382785797 0)),((559440.48582596611231565 5539485.96608389914035797 0, 559430.05711502861231565 5539484.61342276632785797 0, 559441.87690872978419065 5539483.64760245382785797 0, 559440.48582596611231565 5539485.96608389914035797 0)),((559430.99879593681544065 5539497.32625235617160797 0, 559430.05711502861231565 5539484.61342276632785797 0, 559440.48582596611231565 5539485.96608389914035797 0, 559430.99879593681544065 5539497.32625235617160797 0)),((559440.18873734306544065 5539500.42384757101535797 0, 559430.99879593681544065 5539497.32625235617160797 0, 559440.48582596611231565 5539485.96608389914035797 0, 559440.18873734306544065 5539500.42384757101535797 0)),((559430.83259720634669065 5539510.12337149679660797 0, 559430.99879593681544065 5539497.32625235617160797 0, 559440.18873734306544065 5539500.42384757101535797 0, 559430.83259720634669065 5539510.12337149679660797 0)),((559439.89164872001856565 5539514.88161124289035797 0, 559430.83259720634669065 5539510.12337149679660797 0, 559440.18873734306544065 5539500.42384757101535797 0, 559439.89164872001856565 5539514.88161124289035797 0)),((559441.87690872978419065 5539483.64760245382785797 0, 559443.25630326103419065 5539470.04579581320285797 0, 559446.11232132744044065 5539482.49251212179660797 0, 559441.87690872978419065 5539483.64760245382785797 0)),((559446.11232132744044065 5539482.49251212179660797 0, 559443.25630326103419065 5539470.04579581320285797 0, 559453.27732987236231565 5539470.53065909445285797 0, 559446.11232132744044065 5539482.49251212179660797 0)),((559446.11232132744044065 5539482.49251212179660797 0, 559453.27732987236231565 5539470.53065909445285797 0, 559460.78807205986231565 5539482.16267813742160797 0, 559446.11232132744044065 5539482.49251212179660797 0)),((559460.78807205986231565 5539482.16267813742160797 0, 559453.27732987236231565 5539470.53065909445285797 0, 559463.29835648369044065 5539471.01552237570285797 0, 559460.78807205986231565 5539482.16267813742160797 0)),((559440.05949539970606565 5539467.18544425070285797 0, 559443.25630326103419065 5539470.04579581320285797 0, 559429.37440628837794065 5539468.34004630148410797 0, 559440.05949539970606565 5539467.18544425070285797 0)),((559440.05949539970606565 5539467.18544425070285797 0, 559429.37440628837794065 5539468.34004630148410797 0, 559430.51341385673731565 5539467.20106925070285797 0, 559440.05949539970606565 5539467.18544425070285797 0)),((559440.05949539970606565 5539467.18544425070285797 0, 559430.51341385673731565 5539467.20106925070285797 0, 559440.52803177665919065 5539458.28340567648410797 0, 559440.05949539970606565 5539467.18544425070285797 0)),((559440.52803177665919065 5539458.28340567648410797 0, 559430.51341385673731565 5539467.20106925070285797 0, 559430.69905228447169065 5539452.90749014914035797 0, 559440.52803177665919065 5539458.28340567648410797 0)),((559440.52803177665919065 5539458.28340567648410797 0, 559430.69905228447169065 5539452.90749014914035797 0, 559440.99653763603419065 5539449.38136710226535797 0, 559440.52803177665919065 5539458.28340567648410797 0)),((559430.69905228447169065 5539452.90749014914035797 0, 559440.98625321220606565 5539448.98830069601535797 0, 559440.99653763603419065 5539449.38136710226535797 0, 559430.69905228447169065 5539452.90749014914035797 0)))"
        qgs_geom = QgsGeometry()
        qgs_geom = qgs_geom.fromWkt(txt_geom)
        qgs_feature = QgsFeature()
        qgs_feature.setGeometry(qgs_geom)
        ca = ChordalAxis(qgs_feature, GenUtil.ZERO)
        ca.correct_skeleton()
        centre_lines = ca.get_skeleton()

        val0 = True
        self.assertTrue(val0, title)



# Supply path to qgis install location
QgsApplication.setPrefixPath("/usr/bin/qgis", False)

# profile_folder = 'C:\\Users\\berge\\AppData\\Roaming\\QGIS\\QGIS3\\profiles\\test12'
profile_folder = '.'
# Create a reference to the QgsApplication.  Setting the second argument to False disables the GUI.
app = QgsApplication([], False, profile_folder)

# Load providers
app.initQgis()