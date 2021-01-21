from qgis.processing import alg
from qgis.core import QgsFeature, QgsFeatureSink, QgsTopologyPreservingSimplifier, QgsProcessingException,  QgsWkbTypes, QgsGeometry

@alg(name="topological_simplifier", label=alg.tr("Topological simplifier"), group="geosim", group_label=alg.tr("Geo sim"), icon=r"C:\temp\flame.png")
@alg.input(type=alg.SOURCE, name="INPUT", label="Input layer")
@alg.input(type=alg.DISTANCE, name="TOLERANCE", label="Tolerance", default=1.0)
@alg.input(type=alg.BOOL, name="VERBOSE", label="Verbose mode", default=False)
@alg.input(type=alg.SINK, name="OUTPUT", label="Output layer")
@alg.output(type=str, name="NBR_FEATURE", label="Number of features")


def topological_simplifier(instance, parameters, context, feedback, inputs):
    """
    <b>Topological simplifier</b>
    TopoSim is a geospatial simplification tool for lines and polygons. TopoSim implements \
    QGIS's QgsTopologyPreservingSimplifier tool. For line and polygon simplification \
    that tool implements an algorithm similar to the Douglas Peucker algorithm. The implementation \
    preserves the topology within one feature but not between features of the same layer or from \
    different layers. There is also a known bug where the algorithm may create invalid topologies \
    if there are components which are small relative to the tolerance value. In particular, if a \
    small interior hole is very close to an edge, simplification may result in the hole being moved \
    outside the polygon. Toposim will detect these situations where one or more rings (interior \
    parts) fall outside the polygon after being simplified and make the polygon invalid. \
    The algoritm will remove (delete) these ring(s) so the feature remains valid after simplification.

    <b>Usage</b>
    <u>Input</u>: A Line String or Polygon layer
    <u>Tolerance</u>: The tolerance in ground unit used to simplify the line

    For more information: https://github.com/Dan-Eli/GeoSim
    """
    source = instance.parameterAsSource(parameters, "INPUT", context )
    tolerance = instance.parameterAsDouble(parameters,"TOLERANCE", context)
    verbose = instance.parameterAsBool(parameters, "VERBOSE", context)

    if source is None:
        raise QgsProcessingException(instance.invalidSourceError(parameters, "INPUT"))

    (sink, dest_id) = instance.parameterAsSink(parameters, "OUTPUT", context,
                                               source.fields(),
                                               source.wkbType(),
                                               source.sourceCrs()
                                          )

    # Validate input source type
    if source.wkbType() not in [QgsWkbTypes.Polygon, QgsWkbTypes.LineString]:
        raise QgsProcessingException("Can only process: LineString and Polygon type layer")

    # Validate sink
    if sink is None:
        raise QgsProcessingException(instance.invalidSinkError(parameters, "OUTPUT"))

    nbr_feature = 0
    nbr_hole_del = 0
    nbr_feat_inv_correct = 0
    nbr_feat_inv_uncorrect = 0
    total = 100.0 / source.featureCount() if source.featureCount() else 0
    simplifier = QgsTopologyPreservingSimplifier(tolerance)
    features = source.getFeatures()
    for i, feature in enumerate(features):
        geom = feature.geometry()
        s_geom = simplifier.simplify(geom)
        if s_geom.isGeosValid() and s_geom.isSimple():
            # Simplification OK
            pass
        else:
            if s_geom.wkbType() == QgsWkbTypes.Polygon:
                # Process the invalid geometry Polygon
                s_geom_parts = s_geom.asPolygon() # Extract the outer and inner rings in a list
                # Recreate independantPolygon for each part
                s_geom_pols = [QgsGeometry.fromPolygonXY([s_geom_part]) for s_geom_part in s_geom_parts ]
                s_geom_pols.sort(key=polygon_area)  # Sort polygon by ascending area size
                s_geom_outer = s_geom_pols.pop()  # extract the outer ring
                while s_geom_pols:
                    # Extract each inner ring and test if located inside the polygon in construction
                    s_geom_inner = s_geom_pols.pop()
                    if s_geom_inner.within(s_geom_outer):
                        s_geom_outer.addRing(s_geom_inner.asPolygon()[0])  # Add a new inner ring
                    else:
                        nbr_hole_del += 1
                        if verbose:
                            inner_hole = s_geom_inner.asPolygon()
                            inner_hole_xy = inner_hole[0][0]  # Extract the first coordinate of the ring
                            xy = str(inner_hole_xy.x()) + ", " + str(inner_hole_xy.y())
                            feedback.pushInfo("Inner hole deleted: {0}".format(xy))

                if s_geom.isGeosValid() and s_geom.isSimple():
                    # Feature corrected... OK
                    nbr_feat_inv_correct += 1
                else:
                    # Should not happen
                    nbr_feat_inv_uncorrect += 1
            else:
                # Should not happen...
                nbr_feat_inv_uncorrect += 1

        # Add the feature in the sink
        feature.setGeometry(s_geom)
        sink.addFeature(feature, QgsFeatureSink.FastInsert)
        nbr_feature += 1

        if feedback.isCanceled():
            break
        feedback.setProgress(int(i * total))

    # Push some output statistics
    feedback.pushInfo("Number of features simplified: {0}".format(nbr_feature))
    feedback.pushInfo("Number of invalid features corrected: {0}".format(nbr_feat_inv_correct))
    feedback.pushInfo("Number of holes deleted: {0}".format(nbr_hole_del))
    feedback.pushInfo("Number of features left invalid: {0}".format(nbr_feat_inv_uncorrect))

    return {"OUTPUT": dest_id}

def polygon_area(pol):
    return pol.area()