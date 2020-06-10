# Code example available here: https://anitagraser.com/2019/03/03/stand-alone-pyqgis-scripts-with-osgeo4w/

import processing

from qgis.analysis import QgsNativeAlgorithms
from qgis.core import QgsApplication, QgsGeometry, QgsPoint, QgsPointXY, QgsTriangle, QgsFeature, QgsWkbTypes, \
                      QgsSpatialIndex, QgsFeature, QgsRectangle
from processing.core.Processing import Processing
from shapely.geometry import LineString
from datetime import datetime
import random

# Supply path to qgis install location
QgsApplication.setPrefixPath("C:\\OSGEO4~1\\apps\\qgis", True)

#profile_folder = 'C:\\Users\\berge\\AppData\\Roaming\\QGIS\\QGIS3\\profiles\\test12'
profile_folder = '.'
# Create a reference to the QgsApplication.  Setting the second argument to False disables the GUI.
app = QgsApplication([], False, profile_folder)

# Load providers
app.initQgis()

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

loop = 10000
now = datetime.now()

for i in range(loop):
    a = LineString (((0,0),(10,10)))
    c = a.coords[0]


now1 = datetime.now()
print("Elapse: ",now1-now)

index = QgsSpatialIndex()
for i in range(loop):
    geom = QgsGeometry.fromPolylineXY([QgsPointXY(0,0), QgsPointXY(10,10)])
    feat = QgsFeature()
    feat.setGeometry(geom)
    index.addFeature(feat)





0/0
geom = QgsGeometry.fromPolygonXY([[QgsPointXY( 0,  0), QgsPointXY( 5, 10), QgsPointXY(10, 0)],[QgsPointXY( 2,  2), QgsPointXY( 3, 5), QgsPointXY(5, 2)]])


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


def test_spatialindex():
    """This function is testing the spatial index"""

    from shapely.geometry import LineString
    from shapely.strtree import STRtree
    import random, time

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