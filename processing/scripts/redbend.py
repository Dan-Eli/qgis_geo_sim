import os, sys
file = os.path.abspath(__file__)
pathname = os.path.dirname(os.path.abspath(__file__))
sys.path.append(pathname)
from algo_reduce_bend import reduce_bends, is_line_string, is_polygon

from qgis.processing import alg
from qgis.core import QgsFeatureSink, QgsProcessingException, QgsFeatureRequest, QgsWkbTypes

@alg(name="redbend", label=alg.tr("RedBend"), group="geosim", group_label=alg.tr("Geo sim"), icon=r"C:\temp\flame.png")
@alg.input(type=alg.SOURCE, name="INPUT", label="Input layer")
@alg.input(type=alg.DISTANCE, name="DIAMETER", label="Bend diameter", default=1.0)
@alg.input(type=alg.BOOL, name="EXCLUDE_HOLE", label="Exclude holes", default=True)
@alg.input(type=alg.BOOL, name="EXCLUDE_POLYGON", label="Exclude polygons", default=True)
@alg.input(type=alg.SINK, name="OUTPUT", label="RedBend")
@alg.output(type=str, name="NBR_FEATURE", label="Number of features")

def redbend(instance, parameters, context, feedback, inputs):
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

    context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

    source = instance.parameterAsSource(parameters, "INPUT", context )
    exclude_hole = instance.parameterAsBool(parameters, "EXCLUDE_HOLE", context)
    exclude_polygon = instance.parameterAsBool(parameters, "EXCLUDE_POLYGON", context)
    diameter_tol = instance.parameterAsDouble(parameters, "DIAMETER", context)
#    verbose = instance.parameterAsBool(parameters, "VERBOSE", context)

    if source is None:
        raise QgsProcessingException(instance.invalidSourceError(parameters, "INPUT"))

    # Validate input source type
    if is_polygon(source.wkbType()):
        type = QgsWkbTypes.Polygon
    elif is_line_string(source.wkbType()):
        type = QgsWkbTypes.LineString
    else:
        #  Cannot process this feature type
        raise QgsProcessingException("Can only process: LineString or Polygon layers")

    (sink, dest_id) = instance.parameterAsSink(parameters, "OUTPUT", context,
                                               source.fields(),
                                               type,
                                               source.sourceCrs())

    # Validate sink
    if sink is None:
        raise QgsProcessingException(instance.invalidSinkError(parameters, "OUTPUT"))

    features = source.getFeatures()
    qgs_features_in = []
    for qgs_feature_in in features:
        # Load all the QgsFeature
        qgs_features_in.append(qgs_feature_in)

    try:
        # Call the bend reduction
        rb_return = reduce_bends(qgs_features_in, diameter_tol, feedback, exclude_polygon, exclude_hole)
    except Exception:
        import traceback
        traceback.print_exc()

    for qgs_feature_out in rb_return.qgs_features_out:
        sink.addFeature(qgs_feature_out, QgsFeatureSink.FastInsert)

    # Push some output statistics
    feedback.pushInfo("Number of features in: {0}".format(rb_return.in_nbr_features))
    feedback.pushInfo("Number of features out: {0}".format(rb_return.out_nbr_features))
    feedback.pushInfo("Number of iteration needed: {0}".format(rb_return.nbr_pass))
    feedback.pushInfo("Number of bends detected: {0}".format(rb_return.nbr_bends_detected))
    feedback.pushInfo("Number of bends reduced: {0}".format(rb_return.nbr_bends_reduced))
    feedback.pushInfo("Number of deleted polygons: {0}".format(rb_return.nbr_pol_del))
    feedback.pushInfo("Number of deleted polygon holes: {0}".format(rb_return.nbr_hole_del))

    return {"OUTPUT": dest_id, "NBR_FEATURE": rb_return.out_nbr_features}
