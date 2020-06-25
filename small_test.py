# Code example available here: https://anitagraser.com/2019/03/03/stand-alone-pyqgis-scripts-with-osgeo4w/

import processing

from qgis.analysis import QgsNativeAlgorithms
from qgis.core import QgsApplication, QgsGeometry, QgsPoint, QgsPointXY, QgsTriangle, QgsFeature, QgsWkbTypes, \
                      QgsSpatialIndex, QgsFeature, QgsRectangle, QgsLineString
from processing.core.Processing import Processing
from shapely.geometry import LineString
from datetime import datetime
import random
from shapely.strtree import STRtree
import random, time

def use_same_qgspoint():
    """Use the same QgsPoint 2 times"""

    p0 = QgsPointXY(0,0)
    p1 = QgsPointXY(1,1)

    geom0 = QgsGeometry.fromPolylineXY([p0,p1])
    geom1 = QgsGeometry.fromPolylineXY([p0,p1])

    return (geom0,geom1)

def merge_lines():
    """Merge line from a multi line string"""

    geom1 = QgsGeometry.fromPolyline([QgsPoint(0, 0), QgsPoint(10, 10)])
    geom2 = QgsGeometry.fromPolyline([QgsPoint(20, 20), QgsPoint(30, 30)])

    geom_a = geom1.asPolyline()
    geom_b = geom2.asPolyline()

    geoms = [geom_a,geom_b]
    multi_line = QgsGeometry.fromMultiPolylineXY(geoms)

    merge = multi_line.mergeLines()

    for part in merge.constParts():
        print (part)

    return


def test_spatialindex():
    """This function is testing the spatial index of QGIS"""

    lst_lines_shapely = []
    lst_lines_qgis = []
    lst_intersects_shapely = []
    lst_intersects_qgis = []
    # Create the triangles
    for i in range(25000):
        x = random.random() * 10000.
        y = random.random() * 10000.
        coords = [(x, y), (x + 5, y + 5), (x, y + 10), (x, y)]
        lst_lines_shapely.append(LineString(coords))
        geom = QgsGeometry.fromPolylineXY(
            [QgsPointXY(x, y), QgsPointXY(x + 5, y + 5), QgsPointXY(x, y + 10), QgsPointXY(x, y)])
        feat = QgsFeature()
        feat.setGeometry(geom)
        lst_lines_qgis.append(feat)

    # Create the bounding boxes
    for i in range(100000):
        x = random.random() * 10000.
        y = random.random() * 10000.
        coords = [(x, y), (x + 15, y), (x + 15, y + 15), (x, y + 15), (x, y)]
        lst_intersects_shapely.append(LineString(coords))
        rect = QgsRectangle(x, y, x + 15, y + 15)
        lst_intersects_qgis.append(rect)

    # Create shapely STRtree
    tree = STRtree(lst_lines_shapely)

    # Create QGIS index
    index = QgsSpatialIndex()
    index.addFeatures(lst_lines_qgis)

    sec1 = time.time()

    # Find the intersection with STRtree
    str_tree_nbr = 0
    for intersect in lst_intersects_shapely:
        str_tree = tree.query(intersect)
        str_tree_nbr += len(str_tree)

    sec2 = time.time()
    print("Seconds  =", sec2 - sec1)
    print("STRtree number: ", str_tree_nbr)

    # Find the intersection with Qgsspatial Index
    qgis_tree_nbr = 0
    for intersect in lst_intersects_qgis:
        qgis_tree = index.intersects(intersect)
        qgis_tree_nbr += len(qgis_tree)

    sec3 = time.time()
    print("Seconds  =", sec3 - sec2)
    print("STRtree number: ", qgis_tree_nbr)

    return

class GsQgsFeature(QgsFeature):
    # Define a specialization of QgsFeature
    def __init__(self, feature=None):
        super().__init__()
        if feature is not None:
            geom = feature.geometry()
            self.setGeometry(geom)
            null_geom = QgsGeometry.fromWkt('')
            feature.setGeometry(null_geom)

        self._ori_qgs_feature = feature


def super_feature():
    # Test that we can derive a class from from QgsFeature
    feat = GsQgsFeature()
#    feat1 = SuperQgsFeature(QgsFeature())
    geom = QgsGeometry.fromPolylineXY([QgsPointXY(0, 0), QgsPointXY(10, 10), QgsPointXY(20, 20), QgsPointXY(30, 30)])
    feat.setGeometry(geom)
    valid = feat.isValid()
    id = feat.id()
    geom1 = feat.geometry()
    l = geom1.length()

    return

def mid_pnt_line():

    p0 = QgsPointXY(0,0)
    p1 = QgsPointXY(1,1)
    a = QgsGeometry.fromPolylineXY((p0, p1))
    mi = a.interpolate(a.length()/2.)

    return

def test_access_coord():
    """Comparative test speed of accessing Shapely and QGIS coordinate"""

    line1 = QgsLineString([0.,1.], [2.,3.])
    start = line1.startPoint()
    a = line1.asWkt()
    first = line1.startPoint()
    last = line1.endPoint()
    np = line1.numPoints()
    lp = line1.points()
    type1 = line1.wkbType()
    feat = QgsFeature()
    feat.setGeometry(line1)
    line11 = feat.geometry()
    line2 = QgsGeometry.fromPolyline((QgsPoint(0,0), QgsPoint(1,1)))
    gt = line2.get()


    lst_shapely = []
    lst_qgis = []
    nbr_loop = 2000
    nbr_vert = 100
    sec1 = time.time()
    for j in range(nbr_loop):
        coords = ([i+j,i+j] for i in range(nbr_vert))
        line_string = LineString(coords)
        lst_shapely.append(line_string)

    sec2 = time.time()
    print ("Build shapely: ",sec2-sec1)
    for j in range(nbr_loop):
        pxy = (QgsPoint(i+j, i+j) for i in range(nbr_vert))
        geom = QgsGeometry.fromPolyline(pxy)
        lst_qgis.append(geom)

    sec3 = time.time()
    print("Build QGS: ",sec3 - sec2)

    print ("GO***")
    #Test Shapely speed
    l = 0
    sec1 = time.time()
    for line_shapely in lst_shapely:
        dummy = list(line_shapely.coords)
        l += len(dummy)
        dummy= None

    sec2 = time.time()
    print ("Shapely timer: ", sec2-sec1)

    # Test QGIS speed
    l = 0
    for line_qgis in lst_qgis:
        dummy = line_qgis.asPolyline()
        l += len(dummy)
        dummy = None
    sec3 = time.time()
    print("QGIS timer: ", sec3 - sec2)

    # Test QGIS speed
    l = 0
    for line_qgis in lst_qgis:
        line_string = line_qgis.get()
        l = line_string.numPoints()
        dummy = line_string.points()
        l += len(dummy)
        dummy = None
    sec4 = time.time()
    print("QGIS timer: ", sec4 - sec3)


def test_memory_leak():
    """Test that the garbage collector is working well"""

    feats = []
    for i in range(10000):
        pxy = []
        pxy = (QgsPointXY(i,i) for i in range(10000))
        geom = QgsGeometry.fromPolylineXY(pxy)
        feat = QgsFeature()
        feat.setGeometry(geom)
        feats.append(feat)
        geom1 = QgsGeometry.fromWkt('')
        feat.setGeometry(geom1)

    return


def test_create_feature():
    """Create different type of Geometry"""

    # Create a LineString
    geom = QgsGeometry.fromPolylineXY([QgsPointXY(0, 0), QgsPointXY(10, 10), QgsPointXY(20, 20), QgsPointXY(30, 30)])
    feat = QgsFeature()
    feat.setGeometry(geom)
#    feat.setGeometry(QgsGeometry.fromPolylineXY([QgsPointXY(0, 0), QgsPointXY(10, 10), QgsPointXY(20, 20), QgsPointXY(30, 30)]))
    geom = feat.geometry()
    lst_vertice = list(geom.vertices())
    geom_a = geom.constGet()
    first = geom_a[0] # Extract first vertice
    last = geom_a[-1] # Extract last vertice


    # Create a Polygon without a hole
    geom = QgsGeometry.fromPolygonXY([[QgsPointXY(0, 0), QgsPointXY(5, 10), QgsPointXY(10, 0)]])
    feat = QgsFeature()
    feat.setGeometry(geom)
#    feat.setGeometry(QgsGeometry.fromPolygonXY([[QgsPointXY(0, 0), QgsPointXY(5, 10), QgsPointXY(10, 0)]]))
    geom = feat.geometry()
    lst = geom.asPolygon()
    first_last = geom.adjacentVertices(0)
    first_last = geom.adjacentVertices(2)



    # Create a Polygon with a hole
    geom = QgsGeometry.fromPolygonXY([[QgsPointXY(0, 0), QgsPointXY(5, 10), QgsPointXY(10, 0)],
                                      [QgsPointXY(2, 2), QgsPointXY(3, 5), QgsPointXY(5, 2)]])
    feat = QgsFeature()
    feat.setGeometry(geom)

    # extract the information
    geom = feat.geometry()
    polygon = geom.asPolygon()



    # Create a Polygon
    geom = QgsGeometry.fromMultiPolygonXY([[[QgsPointXY(0, 0), QgsPointXY(5, 10), QgsPointXY(10, 0)],
                                      [QgsPointXY(2, 2), QgsPointXY(3, 5), QgsPointXY(5, 2)]]])
    feat = QgsFeature()
    feat.setGeometry(geom)

    #extract the information
    geom = feat.geometry()
    polygon = geom.asPolygon()
    lst_vertice = list(geom.vertices())




    return





# Supply path to qgis install location
QgsApplication.setPrefixPath("C:\\OSGEO4~1\\apps\\qgis", True)

#profile_folder = 'C:\\Users\\berge\\AppData\\Roaming\\QGIS\\QGIS3\\profiles\\test12'
profile_folder = '.'
# Create a reference to the QgsApplication.  Setting the second argument to False disables the GUI.
app = QgsApplication([], False, profile_folder)

# Load providers
app.initQgis()


#-------------------------------------

#Use same points
use_same_qgspoint()

#Merge the lines
merge_lines()

#Find mid point in a line
mid_pnt_line()

# Comparative test speed
test_access_coord()

# Test super class
super_feature()

#test for memory leak with gweometry
test_memory_leak()

# Test creation of different type of feature and geometry
test_create_feature()

# Test the Spatial index
test_spatialindex()




#Create multicollection of TraingleZ

geom1 = QgsGeometry.fromPolylineXY([QgsPointXY(0, 0), QgsPointXY(10, 10)])
feat1 = QgsFeature()
feat1.setGeometry(geom1)
feat1.setId(0)
geom2 = QgsGeometry.fromPolylineXY([QgsPointXY(0, 0), QgsPointXY(10, 10)])
feat2 = QgsFeature()
feat2.setGeometry(geom2)
feat2.setId(10)
index = QgsSpatialIndex()
index.addFeature(feat1)
index.addFeature(feat2)
qgis_tree = index.intersects(QgsRectangle(0,0,10,10))




# Creates some points
point = QgsPoint(5, 10, 0)
feat = QgsFeature()
feat.setGeaometr
point_z = QgsPoint(5, 10, 0)
point_zm = QgsPoint(5, 10, 0, 0)
point = QgsGeometry(point)
point = point.convertToType(QgsWkbTypes.PointGeometry, False)

# Creates a polygonz
points_z = [QgsPoint(0, 0, 0), QgsPoint(5, 10, 0), QgsPoint(10, 0, 0)]
# Create triangles
triangle_a = QgsTriangle(QgsPoint( 0,  0, 0), QgsPoint( 5, 10, 1), QgsPoint(10, 0, 2))
triangle_b = QgsTriangle(QgsPoint(20, 20, 0), QgsPoint(25, 30, 0), QgsPoint(30, 20, 0))

# Extract the coordinate from the triangles
multi_polygon = []
for triangle in [triangle_a,triangle_b]:
    # Extract the vertex from the triangle
    points = [triangle.vertexAt(i) for i in [0,1,2]]
    # Convert the vertex into PointXY
    polygon = [QgsPointXY(point.x(),point.y()) for point in points]
    multi_polygon.append([polygon])

# construc the multipolygon
multi_geom = QgsGeometry.fromMultiPolygonXY(multi_polygon)

for geom in multi_geom.parts():
    print (geom)







geom = QgsGeometry()
print (geom.newPart(triangle_a))
print(geom)
geom = geom.collectGeometry([triangle_a,triangle_b])
geom2 = geom.convertToType(QgsWkbTypes.PolygonGeometry)
print (geom)
print (geom2)
for part in geom2.parts():
    print(part)


#pol = triangle_b.convertToType(QgsWkbTypes.PolygonGeometry, False)
p1 =  triangle_b.vertexAt(0)
pol =  triangle_b.surfaceToPolygon()


feature = QgsFeature()

for part in geom.parts():
    print (part)


ret_code = geom.addPart([triangle_a,triangle_b])
ret_code = geom.addPart(triangle_b)
#ret_code = geom.addPart(part=triangle_a, QgsWkbTypes.GeometryType=QgsWkbTypes.UnknownGeometry)
ret_code = geom.addPart(triangle_b)

for part in geom.parts():
    print (part)









# Stop QGIS appllication
app.exitQgis()
app.exit()


