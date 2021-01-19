import os, sys
file = os.path.abspath(__file__)
pathname = os.path.dirname(os.path.abspath(__file__))
print ("File: ", file)
print ("Path: ", pathname)
sys.path.append(pathname)

from qgis.processing import alg
from qgis.core import QgsFeature, QgsFields, QgsFeatureSink, QgsProcessingException,  QgsWkbTypes, QgsGeometry, \
                      QgsFeatureRequest
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
    if source.wkbType() not in [QgsWkbTypes.MultiPolygon, QgsWkbTypes.MultiPolygonZ]:
        raise QgsProcessingException("Can only process: MultiPolygon and MultiPolygonZ type layer")

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
        # Call the chordal axis
        try:
            ca = ChordalAxis(feature, GenUtil.ZERO)
            if correction:
                ca.correct_skeleton()
            centre_lines = ca.get_skeleton()
        except Exception:
            import traceback
            traceback.print_exc()

        # Load the centre line in the sink
        for line in centre_lines:
            out_feature = QgsFeature()
            geom_feature = QgsGeometry(line.clone())
            out_feature.setGeometry(geom_feature)
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
















