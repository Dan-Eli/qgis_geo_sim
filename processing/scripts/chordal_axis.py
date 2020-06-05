import os, sys
file = os.path.abspath(__file__)
pathname = os.path.dirname(os.path.abspath(__file__))
print ("File: ", file)
print ("Path: ", pathname)
sys.path.append(pathname)

from qgis.processing import alg
from qgis.core import QgsFeature, QgsFields, QgsFeatureSink, QgsProcessingException,  QgsWkbTypes, QgsGeometry, QgsFeatureRequest
from shapely.wkt import dumps, loads
from shapely.geometry import LineString
from lib_geosim import GenUtil, ChordalAxis


@alg(name="chordal_axis", label=alg.tr("Chordal axis"), group="geosim", group_label=alg.tr("Geo sim"), icon=r"C:\temp\flame.png")
@alg.input(type=alg.SOURCE, name="INPUT", label="Input layer")
@alg.input(type=alg.BOOL, name="CORRECTION", label="Correct skeleton", default=True)
@alg.input(type=alg.SINK, name="OUTPUT", label="Output layer")
@alg.output(type=str, name="NBR_FEATURE", label="Number of features")

def chordal_axis(instance, parameters, context, feedback, inputs):
    """
    Given a distance will split a line layer into segments of the distance1
    * coco
    * toto

    <b>This text is bold</b>

    #. Titre
    #. Titre 2
    """

    import os
    print(os.getcwd())
#    context = dataobjects.createContext()
    context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

    source = instance.parameterAsSource(parameters, "INPUT", context )
    correction = instance.parameterAsBool(parameters, "CORRECTION", context)
#    verbose = instance.parameterAsBool(parameters, "VERBOSE", context)

    if source is None:
        raise QgsProcessingException(instance.invalidSourceError(parameters, "INPUT"))

    (sink, dest_id) = instance.parameterAsSink(parameters, "OUTPUT", context,
                                               QgsFields(),
                                               QgsWkbTypes.LineString,
                                               source.sourceCrs()
                                          )

    # Validate input source type
    a = source.wkbType()
    if source.wkbType() not in [QgsWkbTypes.MultiPolygon, QgsWkbTypes.MultiPolygonZ, QgsWkbTypes.Polygon]:
        raise QgsProcessingException("Can only process: MultiPolygon and MultiPolygonZPolygon type layer")

    # Validate sink
    if sink is None:
        raise QgsProcessingException(instance.invalidSinkError(parameters, "OUTPUT"))

    nbr_feature = 0
    nbr_polygon = source.featureCount()
    nbr_triangle = 0
    nbr_centre_line = 0
    total = 100.0 / source.featureCount() if source.featureCount() else 0

    features = source.getFeatures()
    for i, feature in enumerate(features):
        qgis_geom = feature.geometry()
        lst_triangle = []
        for part in qgis_geom.constParts():
            # Transform each part
            shapely_part = qgis_to_shapely(part)
            if len(shapely_part.exterior.coords) != 4:
                0/0
            line = LineString(shapely_part.exterior.coords)
            lst_triangle.append(line)
            nbr_triangle += 1

        # Call the chordal axis
        try:
            ca = ChordalAxis(lst_triangle, GenUtil.ZERO)
            if correction:
                ca.correct_skeleton()
            centre_lines = ca.get_skeleton()
        except Exception:
            import traceback
            traceback.print_exc()

        # Load the centre line in the sink
        for line in centre_lines:
            out_feature = QgsFeature()
            qgis_geom = shapely_to_qgis(line)
            out_feature.setGeometry(qgis_geom)
            sink.addFeature(out_feature, QgsFeatureSink.FastInsert)
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


