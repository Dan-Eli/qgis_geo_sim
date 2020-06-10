import os, sys
file = os.path.abspath(__file__)
pathname = os.path.dirname(os.path.abspath(__file__))
print ("File: ", file)
print ("Path: ", pathname)

sys.path.append(pathname)
from qgis.processing import alg
from qgis.core import QgsFeatureSink, QgsProcessingException,  QgsWkbTypes, QgsGeometry, QgsFeatureRequest
from shapely.wkt import dumps, loads
from algo_sherbend import AlgoSherbend
from lib_geosim import Holder


@alg(name="sherbend", label=alg.tr("Sherbend"), group="geosim", group_label=alg.tr("Geo sim"), icon=r"C:\temp\flame.png")
@alg.input(type=alg.SOURCE, name="INPUT", label="Input layer")
@alg.input(type=alg.DISTANCE, name="DIAMETER", label="Bend diameter", default=1.0)
@alg.input(type=alg.BOOL, name="EXCLUDE_HOLE", label="Exclude holes", default=True)
@alg.input(type=alg.BOOL, name="EXCLUDE_POLYGON", label="Exclude polygons", default=True)
@alg.input(type=alg.SINK, name="OUTPUT", label="Sherbend")
@alg.output(type=str, name="NBR_FEATURE", label="Number of features")

def sherbend(instance, parameters, context, feedback, inputs):
    """
    Sherbend is a geospatial simplification and generalization tool for lines and polygons. \
    Sherbend is an implementation and an improvement of the algorithm described in the paper \
    "Line Generalization Based on Analysis of Shape Characteristics, Zeshen Wang and \
    Jean-Claude Müller, 1998" often known as "Bend Simplify" or "Wang Algorithm". The \
    particularity of this algorithm is that for each line it analyzes its bends (curves) and \
    decides which one to simplify, trying to emulate what a cartographer would do manually \
    to simplify or generalize a line. Sherbend will accept points, lines and polygons as input. \
    Even though points cannot be simplified, they are used for topological relationship \
    validations. Sherbend can accept GeoPackage and Esri Shape file as input/ouput but not a mixed \
    of both.

    <b>Usage</b>

    <u>Input</u> : Any Line string or Polygon layer

    <u>Bend diameter</u>: Theoritical diameter of a bend to remove

    <u>Exclude hole</u>: If you want to exclude holes below the diameter of the bend

    <u>Exclude polygon</u>: If you want to exclude polygon below the diameter of the bend

    <u>Bend diameter</u>: Theoretical diameter of a bend to remove

    <b>Line Simplification versus Line Generalization</b>
    Line Simplification is the process of removing vertices in a line while trying to keep the maximum \
    number of details within the line whereas Line Generalization is the process of removing \
    meaningless (unwanted) details in a line usually for scaling down. The well known Douglas-Peucker \
    algorithm is a very good example of line simplification tool and Sherbend falls more in the \
    category of line generalization tools. Keep in mind thay both algorithms can be complementary because \
    Sherbend will not remove unnecessary vertices in the case of very high densities of vertices. It may \
    be a good idea to use Douglass Peucker before Sherbend in the case of very densed geometries.

    <b>Rule of thumb for the diameter</b>
    Sherbend can be used for line simplifying often in the context of line generalization. The big \
    question will often be what diameter should we use? A good starting point is the cartographic rule of \
    thumb -- the .5mm on the map -- which says that the minimumm distance between two lines should be \
    greater than 0.5mm on a paper map. So to simplify (generalize) a line for representation at a scale of \
    1:50 000 for example a diameter of 25m should be a good starting point...Rule of thumb for the diameter \
    Sherbend can be used for line simplifying often in the context of line generalization. The big question will \
    often be what diameter should we use? A good starting point is the cartographic rule of thumb -- the .5mm \
    on the map -- which says that the minimumm distance between two lines should be greater than 0.5mm on a paper \
    map. So to simplify (generalize) a line for representation at a scale of 1:50 000 for example a diameter of \
    25m should be a good starting point...

    for more information: https://github.com/Dan-Eli/GeoSim

    """

#    context = dataobjects.createContext()
    context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

    source = instance.parameterAsSource(parameters, "INPUT", context )
    exclude_hole = instance.parameterAsBool(parameters, "EXCLUDE_HOLE", context)
    exclude_polygon = instance.parameterAsBool(parameters, "EXCLUDE_POLYGON", context)
    diameter = instance.parameterAsDouble(parameters, "DIAMETER", context)
#    verbose = instance.parameterAsBool(parameters, "VERBOSE", context)

    if source is None:
        raise QgsProcessingException(instance.invalidSourceError(parameters, "INPUT"))

    # Validate input source type
    if source.wkbType() not in [QgsWkbTypes.LineString, QgsWkbTypes.Polygon]:
        # Can only process LineString or Polygon
        raise QgsProcessingException("Can only process: LineString or Polygon type layer")

    (sink, dest_id) = instance.parameterAsSink(parameters, "OUTPUT", context,
                                               source.fields(),
                                               source.wkbType(),
                                               source.sourceCrs()
                                              )

    layer_name = 'layer_name'

    dlayer_dict = {layer_name:diameter}
    command = Holder(exclude_hole=exclude_hole, exclude_polygon=exclude_polygon, dlayer_dict=dlayer_dict)
    geo_content = Holder( in_nbr_points=0, in_nbr_line_strings=0, in_nbr_polygons=0, in_nbr_holes=0,
                                  out_nbr_points=0, out_nbr_line_strings=0, out_nbr_polygons=0, out_nbr_holes=0,
                                  nbr_del_polygons=0, nbr_del_holes=0, nbr_bends_simplified=0)
    # Validate sink
    if sink is None:
        raise QgsProcessingException(instance.invalidSinkError(parameters, "OUTPUT"))

    total = 100.0 / source.featureCount() if source.featureCount() else 0

    features = source.getFeatures()
    for i, qgis_feat in enumerate(features):
        qgis_geom = qgis_feat.geometry()
        if qgis_geom.wkbType() == QgsWkbTypes.LineString:
            geo_content.in_nbr_line_strings += 1
        else:
            geo_content.in_nbr_polygons += 1
        qgis_geom = qgis_feat.geometry()
        shapely_feat = qgis_to_shapely(qgis_geom)
        shapely_feat.sb_layer_name = layer_name
        shapely_feat.sb_properties = None
        geo_content.in_features = [shapely_feat]
        geo_content.out_features = []

        try:
            sherbend = AlgoSherbend(command, geo_content)
            sherbend.process()
        except Exception:
            import traceback
            traceback.print_exc()

        for shapely_feat in geo_content.out_features:
            qgis_geom = shapely_to_qgis(shapely_feat)
            if qgis_geom.wkbType() == QgsWkbTypes.LineString:
                geo_content.out_nbr_line_strings += 1
            else:
                geo_content.out_nbr_polygons += 1
            qgis_feat.setGeometry(qgis_geom)
            sink.addFeature(qgis_feat, QgsFeatureSink.FastInsert)

        if feedback.isCanceled():
            break

        feedback.setProgress(int(i * total))


    # Push some output statistics
    feedback.pushInfo("Number of in line strings: {0}".format(geo_content.in_nbr_line_strings))
    feedback.pushInfo("Number of in polygons: {0}".format(geo_content.in_nbr_polygons))
    feedback.pushInfo("Number of out line strings: {0}".format(geo_content.out_nbr_line_strings))
    feedback.pushInfo("Number of out polygons: {0}".format(geo_content.out_nbr_polygons))
    feedback.pushInfo("Number of deleted polygons: {0}".format(geo_content.nbr_del_polygons))
    feedback.pushInfo("Number of deleted polygon holes: {0}".format(geo_content.nbr_del_holes))
    feedback.pushInfo("Number of bends simplified: {0}".format(geo_content.nbr_bends_simplified))

    return {"OUTPUT": dest_id, "NBR_FEATURE": geo_content.out_nbr_polygons}


def shapely_to_qgis(shapely_geom):

    wkt = dumps(shapely_geom)
    qgis_geom = QgsGeometry.fromWkt(wkt)

    return qgis_geom


def qgis_to_shapely(qgis_geom):

    wkt = qgis_geom.asWkt()
    shapely_geom = loads(wkt)

    return shapely_geom


