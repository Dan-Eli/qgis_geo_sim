import os, sys
file = os.path.abspath(__file__)
pathname = os.path.dirname(os.path.abspath(__file__))
sys.path.append(pathname)

from qgis.processing import alg
from qgis.core import QgsFeature, QgsFields, QgsFeatureSink, QgsProcessingException,  QgsWkbTypes, QgsGeometry, \
                      QgsFeatureRequest, QgsPointXY
from shapely.wkt import dumps, loads
from shapely.geometry import LineString
from lib_geosim import GenUtil, ChordalAxis


@alg(name="chordal_axis", label=alg.tr("Chordal axis"), group="geosim", group_label=alg.tr("Geo sim"), icon=r"C:\temp\flame.png")
@alg.input(type=alg.SOURCE, name="INPUT", label="Input layer")
@alg.input(type=alg.BOOL, name="CORRECTION", label="Correct skeleton", default=True)
@alg.input(type=alg.SINK, name="OUTPUT", label="Chordal axis")
@alg.output(type=str, name="NBR_FEATURE", label="Number of features")

def chordal_axis(instance, parameters, context, feedback, inputs):
    """
    <b>Chordal Axis</b>
    ChordalAxis is a geospatial tool that takes triangles, usually the result of a constraint \
    Delauny trianglulation and creates a skeleton (the center line). ChordalAxis is an improvement \
    of the algorithm based of the paper "Rectification of the Chordal Axis Transform and a \
    New Criterion for Shape Decomposition", Lakshman Prasad, 2005".

    <b>Medial Axis Versus Chordal Axis</b>
    The skeleton (center line) is a linear feature representation of a polygonized feature. In \
    computational geometry, it is known as the medial axis and many algorithms are approximating \
    it very well. A major issue with those algorithms is the possible instability for very irregular \
    complex polygons such as dense river or road network polygons. (Figure 4). The Chordal Axis has \
    shown excellent stability in very very polygons while extracting a very representative skeleton.

    <b>Usage</b>
    <u>Input</u>: A Multipolygon layer where each triangle composing a polygon is part of the same multipolygon. \
    For example the  output of Tessellate processing script.

    <u>Correct skeleton</u>:  Correct the skeleton for small centre line, T junction and X junction. Usefull in the case \
    of long any narrow polygon (ex.: polygonized road network)

    For more information: https://github.com/Dan-Eli/GeoSim
    """

    # In case an invalid geometry is in the file
    context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

    # read the parameters
    source = instance.parameterAsSource(parameters, "INPUT", context )
    correction = instance.parameterAsBool(parameters, "CORRECTION", context)
#    verbose = instance.parameterAsBool(parameters, "VERBOSE", context)

    if source is None:
        raise QgsProcessingException(instance.invalidSourceError(parameters, "INPUT"))

    (sink, dest_id) = instance.parameterAsSink(parameters, "OUTPUT", context,
                                               QgsFields(),
                                               QgsWkbTypes.LineString,
                                               source.sourceCrs() )

    if source.wkbType() not in [QgsWkbTypes.MultiPolygon, QgsWkbTypes.MultiPolygonZ]:
        raise QgsProcessingException("Can only process: MultiPolygon and MultiPolygonZ type layer")

    # Validate sink
    if sink is None:
        raise QgsProcessingException(instance.invalidSinkError(parameters, "OUTPUT"))

    nbr_polygon = source.featureCount()
    nbr_triangle = 0
    nbr_centre_line = 0
    total = 100.0 / source.featureCount() if source.featureCount() else 0

    features = source.getFeatures()
    for i, multi_feature in enumerate(features):
#        qgis_geom = feature.geometry()
#        lst_triangle = []
#
#        for part in qgis_geom.constParts():
#            # Transform each part into polygon
#            part = to_polygon(part)
#            shapely_part = qgis_to_shapely(part)
#            if len(shapely_part.exterior.coords) != 4:
#                raise QgsProcessingException("A triangle polygon must have exactly four vertice")
#            line = LineString(shapely_part.exterior.coords)
#            lst_triangle.append(line)
#            nbr_triangle += 1

        # Call the chordal axis
        try:
            ca = ChordalAxis(multi_feature, GenUtil.ZERO)
            if correction:
                ca.correct_skeleton()
            centre_lines = ca.get_skeleton()
        except Exception:
            import traceback
            traceback.print_exc()

        # Load the centre line in the sink
        for qgis_line in centre_lines:
#            out_feature = QgsFeature()
#            qgis_geom = shapely_to_qgis(line)
#            out_feature.setGeometry(qgis_geom)
            try:
                sink.addFeature(qgis_line, QgsFeatureSink.FastInsert)
            except Exception:
                import traceback
                traceback.print_exc()
            nbr_centre_line + 1

        if feedback.isCanceled():
            break

        feedback.setProgress(int(i * total))

    # Push some output statistics
    feedback.pushInfo("Number of polygons: {0}".format(nbr_polygon))
    feedback.pushInfo("Number of triangles: {0}".format(nbr_triangle))
    feedback.pushInfo("Number of centre_lines: {0}".format(nbr_centre_line))

    return {"OUTPUT": dest_id, "NBR_FEATURE": nbr_centre_line}


def shapely_to_qgis(shapely_geom):

    wkt = dumps(shapely_geom)
    qgis_geom = QgsGeometry.fromWkt(wkt)

    return qgis_geom


def qgis_to_shapely(qgis_geom):

    wkt = qgis_geom.asWkt()
    shapely_geom = loads(wkt)

    return shapely_geom

def to_polygon(qgis_geom):
    # Extract the coordinate from the triangles

    if qgis_geom.wkbType() in (QgsWkbTypes.Triangle, QgsWkbTypes.TriangleM,
                             QgsWkbTypes.TriangleZ, QgsWkbTypes.TriangleZM):
        # Extract the vertex from the triangle
        points = [qgis_geom.vertexAt(i) for i in [0, 1, 2]]
        # Convert the vertex into PointXY
        points_xy= [QgsPointXY(point.x(), point.y()) for point in points]
        # Reconstrucu the polygon
        qgis_geom = QgsGeometry.fromPolygonXY([points_xy])

    elif qgis_geom.wkbType() == QgsWkbTypes.Polygon:
        pass
    else:
        raise QgsProcessingException("Cannot process WKB Type: {0}".format(qgis_geom.wkbType()))

    return qgis_geom






