import os, sys
file = os.path.abspath(__file__)
pathname = os.path.dirname(os.path.abspath(__file__))
print ("File: ", file)
print ("Path: ", pathname)

sys.path.append(pathname)
from qgis.processing import alg
from qgis.core import QgsFeatureSink, QgsProcessingException, QgsFeatureRequest

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

#    context = dataobjects.createContext()
    context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

    source = instance.parameterAsSource(parameters, "INPUT", context )
    exclude_hole = instance.parameterAsBool(parameters, "EXCLUDE_HOLE", context)
    exclude_polygon = instance.parameterAsBool(parameters, "EXCLUDE_POLYGON", context)
    diameter_tol = instance.parameterAsDouble(parameters, "DIAMETER", context)
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

#    dlayer_dict = {layer_name:diameter}
#    command = Holder(exclude_hole=exclude_hole, exclude_polygon=exclude_polygon, dlayer_dict=dlayer_dict)
#    geo_content = Holder( in_nbr_points=0, in_nbr_line_strings=0, in_nbr_polygons=0, in_nbr_holes=0,
#                                  out_nbr_points=0, out_nbr_line_strings=0, out_nbr_polygons=0, out_nbr_holes=0,
#                                  nbr_del_polygons=0, nbr_del_holes=0, nbr_bends_simplified=0)
    # Validate sink
    if sink is None:
        raise QgsProcessingException(instance.invalidSinkError(parameters, "OUTPUT"))

    total = 100.0 / source.featureCount() if source.featureCount() else 0

    features = source.getFeatures()
    qgs_features_in = []
    for qgs_feature_in in features:
        qgs_features_in.append(qgs_feature_in)
#        qgis_geom = qgis_feat.geometry()
#        if qgis_geom.wkbType() == QgsWkbTypes.LineString:
#            geo_content.in_nbr_line_strings += 1
#        else:
#            geo_content.in_nbr_polygons += 1
#        qgis_geom = qgis_feat.geometry()
#        shapely_feat = qgis_to_shapely(qgis_geom)
#        shapely_feat.sb_layer_name = layer_name
#        shapely_feat.sb_properties = None
#        geo_content.in_features = [shapely_feat]
#        geo_content.out_features = []

    try:
        rb_return = reduce_bends(qgs_features_in, diameter_tol)
#            sherbend = AlgoSherbend(command, geo_content)
#            sherbend.process()
    except Exception:
        import traceback
        traceback.print_exc()

#        for shapely_feat in geo_content.out_features:
#            qgis_geom = shapely_to_qgis(shapely_feat)
#            if qgis_geom.wkbType() == QgsWkbTypes.LineString:
#                geo_content.out_nbr_line_strings += 1
#            else:
#                geo_content.out_nbr_polygons += 1
#            qgis_feat.setGeometry(qgis_geom)
    for qgs_feature_out in rb_return.qgs_features_out:
        sink.addFeature(qgs_feature_out, QgsFeatureSink.FastInsert)

#    if feedback.isCanceled():
#        break

#        feedback.setProgress(int(i * total))

    # Push some output statistics
    feedback.pushInfo("Number of features in: {0}".format(rb_return.in_nbr_features))
    feedback.pushInfo("Number of features out: {0}".format(rb_return.out_nbr_features))
    feedback.pushInfo("Number of bends detected: {0}".format(rb_return.nbr_bends_detected))
    feedback.pushInfo("Number of bends reduced: {0}".format(rb_return.nbr_bends_reduced))
#    feedback.pushInfo("Number of deleted polygons: {0}".format(geo_content.nbr_del_polygons))
#    feedback.pushInfo("Number of deleted polygon holes: {0}".format(geo_content.nbr_del_holes))

#    return {"OUTPUT": dest_id, "NBR_FEATURE": geo_content.out_nbr_polygons}
    return {"OUTPUT": dest_id, "NBR_FEATURE": rb_return.out_nbr_features}










# Remaining modifications:
#  - attribute for isClosed instead off call each time isClosed()
#  - Manage case where a closed line string is of the same orientation except the start/end which is different
#  - Have a dynamic epsilon for some cases
#  - Add some comprehension list
#  - Test performance with the profiler
#  - Manage delete Polygon and/or holes below diameter tolerance
from abc import ABC, abstractmethod
import math
from qgis.core import QgsLineString, QgsPolygon, QgsWkbTypes, QgsSpatialIndex,\
                      QgsGeometry, QgsGeometryUtils, QgsVertexId

# Define global constant
GeoSim_EPSILON = 1.0E-9
GeoSim_CW = 0
GeoSim_ACW = -1


class RbFeature(ABC):

    _id_counter = 0

    def __init__(self, qgs_feature):
        self.qgs_feature = qgs_feature
        self.id = RbFeature._id_counter
        RbFeature._id_counter += 1
        abs_geom = qgs_feature.geometry().constGet()
        self.qgs_geom = QgsGeometry(abs_geom.clone())
        self.qgs_feature.clearGeometry()

    @abstractmethod
    def get_rb_geom(self):
        pass

    @abstractmethod
    def get_qgs_feature(self):
        pass


class RbPolygon(RbFeature):

    def __init__(self, qgs_feature):
        super().__init__(qgs_feature)

        if self.qgs_geom.wkbType() != QgsWkbTypes.Polygon:
            # Transform geometry to Polygon
            self.qgs_geom = self.qgs_geom.coerceToType(QgsWkbTypes.Polygon)

        # Transform geometry into a list a LineString first ring being outer ring
        self.qgs_geom = self.qgs_geom.coerceToType(QgsWkbTypes.LineString)
        self.rb_geom = [RbGeom(qgs_geom) for qgs_geom in self.qgs_geom]
        self.qgs_geom = None

    def get_rb_geom(self):
        return self.rb_geom

    def get_qgs_feature(self):
        qgs_pol = QgsPolygon()
        qgs_pol.setExteriorRing(self.rb_geom[0].qgs_geom.constGet().clone())
        for rb_geom in self.rb_geom[1:]:
            qgs_pol.addInteriorRing(rb_geom.qgs_geom.constGet().clone())
        self.qgs_feature.setGeometry(qgs_pol)

        return self.qgs_feature


class RbLineString(RbFeature):

    def __init__(self, qgs_feature):
        super().__init__(qgs_feature)

        if self.qgs_geom.wkbType() != QgsWkbTypes.LineString:
            # Transform geometry to Point
            self.qgs_geom = self.qgs_geom.coerceToType(QgsWkbTypes.LineString)

        self.rb_geom = [RbGeom(self.qgs_geom)]
        self.qgs_geom = None

    def get_rb_geom(self):

        return self.rb_geom

    def get_qgs_feature(self):
        qgs_geom = QgsGeometry(self.rb_geom[0].qgs_geom.constGet().clone())
        self.qgs_feature.setGeometry(qgs_geom)
        return self.qgs_feature


class RbPoint(RbFeature):

    def __init__(self, qgs_feature):
        super().__init__(qgs_feature)

        if self.qgs_geom.wkbType() != QgsWkbTypes.Point:
            # Transform geometry to QgsPoint
            self.qgs_geom = self.qgs_geom.coerceToType(QgsWkbTypes.Point)

        self.rb_geom = [RbGeom(self.qgs_geom)]
        self.qgs_geom = None

    def get_rb_geom(self):
        return self.rb_geom

    def get_qgs_feature(self):
        qgs_geom = QgsGeometry(self.rb_geom[0].qgs_geom.constGet().clone())
        self.qgs_feature.setGeometry(qgs_geom)
        return self.qgs_feature


class RbCollection(object):

    def __init__(self):

        self._spatial_index = QgsSpatialIndex()
        self.dict_rb_features = {}

        return

    def add_features(self, rb_geoms):

        for rb_geom in rb_geoms:
            self.dict_rb_features[rb_geom.id] = rb_geom.qgs_geom
            self._spatial_index.addFeature(rb_geom.id, rb_geom.qgs_rectangle)

        return

    def get_features(self, qgs_rectangle, drop_ids=[]):

        keys = self._spatial_index.intersects(qgs_rectangle)
        # Transform the keys into rb_features
        rb_features = [self.dict_rb_features[key] for key in keys if key not in drop_ids]

        return rb_features


class RbGeom:

    _id_counter = 0

    def __init__(self, qgs_geom):

        self.id = RbGeom._id_counter
        self.qgs_geom = QgsGeometry(qgs_geom)
        self.qgs_rectangle = self.qgs_geom.boundingBox()
        self.bends = []
        self.id = RbGeom._id_counter
        RbGeom._id_counter += 1
        self.is_simplest = False
        self.nbr_bend_reduced = 0
        self.nbr_bend_detected = None

    def edit_closed_line(self, diameter_tol):

        if self.qgs_geom.wkbType() == QgsWkbTypes.LineString:
            bend_area_ok = None
            qgs_line_string = self.qgs_geom.constGet()
            if qgs_line_string.isClosed():
                for bend in self.bends:
                    if bend.bend_area >= calculate_min_adj_area(diameter_tol):
                        bend_area_ok = bend
                        if bend.j - bend.i >= 4:
                            bend_area_ok = bend
                            break
                if bend_area_ok is not None:
                    # There is bend candidate for a rotation
                    qgs_points = qgs_line_string.points()
                    new_start_end = (bend_area_ok.j+bend_area_ok.i) // 2
                    new_qgs_points = qgs_points[new_start_end:] + qgs_points[1:new_start_end+1]
                    self.qgs_geom = QgsGeometry(QgsLineString(new_qgs_points))
                    detect_bends(self)


class Bend:

    def __init__(self, i, j, qgs_points):

        self.i = i
        self.j = j
        self.qgs_geom_bend = QgsPolygon(QgsLineString(qgs_points[i:j+1]))
        self.bend_area = self.qgs_geom_bend.area()
        self.bend_perimeter = self.qgs_geom_bend.perimeter()
        self.compact_index = 4 * self.bend_area * math.pi / self.bend_perimeter**2
        self.adj_area = self.bend_area * (.75/self.compact_index)
        self.to_reduce = False
        self.qgs_geom_new_subline = None
        self.qgs_geom_new_subline_trimmed = None

    def reduce(self, rb_geom):

        for ind in range(self.j-1, self.i, -1):  # Process the vertice in reverse order
            rb_geom.qgs_geom.deleteVertex(ind)  # Delete vertex to reduce line

        # Update statistics
        rb_geom.nbr_bend_reduced += 1

    def get_new_subline(self, rb_geom):

        if self.qgs_geom_new_subline is None:
            # First pass calculate the value
            qgs_pnt_i = rb_geom.qgs_geom.vertexAt(self.i)
            qgs_pnt_j = rb_geom.qgs_geom.vertexAt(self.j)
            qgs_ls_new_subline = QgsLineString([qgs_pnt_i, qgs_pnt_j])
            self.qgs_geom_new_subline = QgsGeometry(qgs_ls_new_subline)

        return self.qgs_geom_new_subline

    def get_new_subline_trimmed(self, rb_geom):

        if self.qgs_geom_new_subline_trimmed is None:
            qgs_ls_new_line = self.get_new_subline(rb_geom).constGet()
            qgs_pnt_i_trimmed = qgs_ls_new_line.interpolatePoint(GeoSim_EPSILON)
            qgs_pnt_j_trimmed = qgs_ls_new_line.interpolatePoint(qgs_ls_new_line.length() - GeoSim_EPSILON)
            qgs_ls_new_subline_trimmed = QgsLineString([qgs_pnt_i_trimmed, qgs_pnt_j_trimmed])
            self.qgs_geom_new_subline_trimmed = QgsGeometry(qgs_ls_new_subline_trimmed)

        return self.qgs_geom_new_subline_trimmed


class RbResults:
    
    def __init__(self):

        in_nbr_features = None
        out_nbr_features = None
        nbr_bend_reduced = None
        qgs_features_out = None
    
    

def dummy():
   """...

    Parameters
    ----------
    aaa : str
        blblbla

    Returns
    -------
    list
       blablabla
    """


def _create_rb_feature(qgs_features):
    """Break the QgsPolygon into RbLineString (inner and outer ring)

    If the feature is a QgsFeature the feature is transformed into a RbLineString

     Parameters
     ----------
     qgs_features : List of QgsFeatures
         List of QgsFeature to simplify

     Returns
     -------
     None
     """

    rb_features = []

    for qgs_feature in qgs_features:
        geom = qgs_feature.geometry()  # extract the Geometry

        if geom.wkbType() in [QgsWkbTypes.Polygon,
                              QgsWkbTypes.Polygon25D,
                              QgsWkbTypes.PolygonZ,
                              QgsWkbTypes.PolygonM,
                              QgsWkbTypes.PolygonZM]:
            rb_features.append(RbPolygon(qgs_feature))

        elif geom.wkbType() in [QgsWkbTypes.LineString,
                                QgsWkbTypes.LineString25D,
                                QgsWkbTypes.LineStringZ,
                                QgsWkbTypes.LineStringM,
                                QgsWkbTypes.LineStringZM]:
            rb_features.append(RbLineString(qgs_feature))

        elif geom.wkbType() in [QgsWkbTypes.Point,
                                QgsWkbTypes.Point25D,
                                QgsWkbTypes.PointM,
                                QgsWkbTypes.PointZ,
                                QgsWkbTypes.PointZM]:
            rb_features.append(RbPoint(qgs_feature))
        else:
            print("Raise exception unsupported GeometryType...")

    return rb_features


def pseudo_reduce(rb_geom):

    # Delete almost duplicate point
    qgs_line_string = rb_geom.qgs_geom.constGet()
    qgs_line_string.removeDuplicateNodes(GeoSim_EPSILON)

    # Delete co-linear points on the line angle of near 0 or 180 degrees
    qgs_points = qgs_line_string.points()
    num_points = len(qgs_points)
    i = 1
    while i <= num_points-2 and num_points >= 3:
        p0_x, p0_y = qgs_points[i-1].x(), qgs_points[i-1].y()
        p1_x, p1_y = qgs_points[i].x(), qgs_points[i].y()
        p2_x, p2_y = qgs_points[i+1].x(), qgs_points[i+1].y()
        angle = QgsGeometryUtils.angleBetweenThreePoints(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y)
        if abs(angle-math.pi) <= GeoSim_EPSILON or abs(angle) <= GeoSim_EPSILON:
            # Co-linear point or flat angle delete the current point
            del qgs_points[i]
            vertex_id = QgsVertexId(_part=0, _ring=0, _vertex=i, _type=QgsVertexId.SegmentVertex)
            qgs_line_string.deleteVertex(vertex_id)
            num_points -= 1
        else:
            i += 1

    return


def detect_bends(rb_geom):

    # Delete almost duplicate point
    qgs_line_string = rb_geom.qgs_geom.constGet()
    qgs_points = qgs_line_string.points()
    num_points = len(qgs_points)
    i = 1
    angles = []
    # Extract the angles at each vertice except start and end vertice
    while i <= num_points - 2 and num_points >= 3:
        p0_x, p0_y = qgs_points[i - 1].x(), qgs_points[i - 1].y()
        p1_x, p1_y = qgs_points[i].x(), qgs_points[i].y()
        p2_x, p2_y = qgs_points[i + 1].x(), qgs_points[i + 1].y()
        angles.append(QgsGeometryUtils.angleBetweenThreePoints(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y))
        i += 1

    # Modify the angle to binary value: clockwise or anti clockwise
    orientation = [GeoSim_CW if angle >= math.pi else GeoSim_ACW for angle in angles]
    rb_geom.bends = []  # Reset the list of bends

    if len(orientation) >= 1:
        start = 0
        direction = orientation[0]
        i = 1
        while i < len(orientation):
            if direction == orientation[i]:
                pass   # Nothing to do
            else:
                # Change of direction.  Create a bend
                rb_geom.bends.append(Bend(start, i + 1, qgs_points))
                direction = orientation[i]
                start = i
            i += 1  # Go look for the next vertice

        rb_geom.bends.append(Bend(start, i + 1, qgs_points))  # Create the last bend

        if qgs_line_string.isClosed() and \
           i+2-start == num_points:
            # A closed lines with only one bend cannot be reduced
            rb_geom.bends = []  # Reset the bends

    if len(rb_geom.bends) == 0:
        # There is no bend so nothing to simplify
        rb_geom.is_simplest = True

    if rb_geom.nbr_bend_detected is None:
        # For the first pass set the number of bend detected (for statistics purpose)
        rb_geom.nbr_bend_detected = len(rb_geom.bends)

    return


def _rb_geom_to_reduce(last_nbr_bend_reduced, rb_geoms):

    new_nbr_bend_reduced = sum(rb_geom.nbr_bend_reduced for rb_geom in rb_geoms)
    if new_nbr_bend_reduced == last_nbr_bend_reduced:
        bend_to_reduce = False
    else:
        bend_to_reduce = True

    return bend_to_reduce


def remove_rb_geoms_done(rb_geoms, rb_geoms_done):

    to_move = []
    for i, rb_geom in enumerate(rb_geoms):
        if rb_geom.qgs_geom.wkbType() == QgsWkbTypes.Point or \
           rb_geom.is_simplest:
            to_move.append(i)

    # Move from rb_geoms to rb_geoms_done
    for i in to_move[::-1]:  # Process the list in reverse order
        rb_geoms_done.append(rb_geoms[i])
        del rb_geoms[i]

    return rb_geoms


def calculate_min_adj_area(diameter_tol):

    min_adj_area = .75 * math.pi * (diameter_tol/2.)**2

    return min_adj_area


def flag_bend_to_reduce(rb_geom, diameter_tol):

    # Minimum adjusted area used to find bend to reduce
    min_adj_area = .75 * math.pi * (diameter_tol/2.)**2

    if rb_geom.qgs_geom.constGet().isClosed and len(rb_geom.bends) >= 3:
        # The closed line has been rotated and the start/end point lie on a bend that do not need to be reduced
        del rb_geom.bends[0]   # Remove the first bend
        del rb_geom.bends[-1]  # Remove the last bend

    lst_bends = []
    for i, bend in enumerate(rb_geom.bends):
        lst_bends.append((bend.adj_area, i))

    # Sort from smallest to biggest bend area
    lst_bends.sort(key=lambda item: item[0])
    start = 0
    end = len(rb_geom.bends)-1
    at_least_one = False

    for (adj_area, i) in lst_bends:
        if adj_area <= min_adj_area:
            if len(lst_bends) == 1:
                rb_geom.bends[i].to_reduce = True  # Only one bend process it...
            else:
                if i == start:
                    if rb_geom.bends[i+1].to_reduce:
                        pass  # Cannot reduce two bend adjacent
                    else:
                        rb_geom.bends[i].to_reduce = True
                        at_least_one = True
                elif i == end:
                    if rb_geom.bends[i-1].to_reduce:
                        pass  # Cannot reduce two bend adjacent
                    else:
                        rb_geom.bends[i].to_reduce = True
                        at_least_one = True
                elif rb_geom.bends[i-1].to_reduce or rb_geom.bends[i+1].to_reduce:
                    pass  # Cannot reduce two bend adjacent
                else:
                    rb_geom.bends[i].to_reduce = True
        else:
            # Over minimum adjusted area
            break

        if not at_least_one:
            rb_geom.is_simplest = True


def validate_spatial_constraints(ind, rb_geom, rb_collection):

    check_constraints = True
    bend = rb_geom.bends[ind]
    qgs_geom_new_subline = bend.get_new_subline(rb_geom)
    qgs_geom_new_subline_trimmed = bend.get_new_subline_trimmed(rb_geom)
    qgs_geom_engine_new_subline = QgsGeometry.createGeometryEngine(qgs_geom_new_subline.constGet())
    qgs_geom_engine_new_subline.prepareGeometry()
    qgs_geom_engine_new_subline_trimmed = QgsGeometry.createGeometryEngine(qgs_geom_new_subline_trimmed.constGet())
    qgs_geom_line_string = rb_geom.qgs_geom

    # First: check if the bend reduce line string is an OGC simple line
    # We test with a tiny smaller line to ease the testing and false positive error
    if qgs_geom_engine_new_subline_trimmed.intersects(qgs_geom_line_string.constGet()):
        # The new sub line intersect the line itself. The result would create a non OGC simple line
        check_constraints = False

    # Second: check that the new line does not intersect any other line or points
    if check_constraints:
        qgs_geom_engine_bend_area = QgsGeometry.createGeometryEngine(bend.qgs_geom_bend)
        qgs_rectangle = bend.qgs_geom_bend.boundingBox()
        qgs_geom_potentials = rb_collection.get_features(qgs_rectangle, [rb_geom.id])
        for qgs_geom_potential in qgs_geom_potentials:
            if qgs_geom_engine_new_subline.intersects(qgs_geom_potential.constGet()):
                # The bend area intersects with a point
                check_constraints = False
                break

    # Third: check that inside the bend to reduce there is no feature completely inside it.  This would cause a
    # sidedness or relative position error
    if check_constraints:
        for qgs_geom_potential in qgs_geom_potentials:
            if qgs_geom_engine_bend_area.contains(qgs_geom_potential.constGet()):
                # A feature is totaly located inside
                check_constraints = False
                break

    return check_constraints


def process_bends(rb_geom, rb_collection):

    for ind in reversed(range(len(rb_geom.bends))):
        bend = rb_geom.bends[ind]
        if bend.to_reduce:
            # Check spatial constraints
            spatial_constraints = validate_spatial_constraints(ind, rb_geom, rb_collection)
            if spatial_constraints:
                bend.reduce(rb_geom)


def _manage_reduce_bend(rb_geoms, rb_collection, diameter_tol):

    rb_geoms_done = []
    nbr_pass = 0
    bend_to_reduce = True
    while bend_to_reduce:

        rb_geoms_done = remove_rb_geoms_done(rb_geoms, rb_geoms_done)  # Remove item to accelerate process
        last_nbr_bend_reduced = sum(rb_geom.nbr_bend_reduced for rb_geom in rb_geoms)
        for rb_geom in rb_geoms:
            pseudo_reduce(rb_geom)  # Remove co-linear and almost touching points
            detect_bends(rb_geom)
            if nbr_pass == 0:
                # Edit the start/end point for closed QgsLineString
                rb_geom.edit_closed_line(diameter_tol)
            flag_bend_to_reduce(rb_geom, diameter_tol)
            process_bends(rb_geom, rb_collection)
        # Check if all bend are processed
        bend_to_reduce = _rb_geom_to_reduce(last_nbr_bend_reduced, rb_geoms)
        nbr_pass += 1

    # Reset the rb_geoms list
    rb_geoms += rb_geoms_done

    return


def reduce_bends(qgs_in_features, diameter_tol):
    """Main routine for bend reduction

     Parameters
     ----------
     qgs_in_features : List of QgsFeatures
         List of QgsFeatures to reduce
    diameter_tol : Float
         Tolerance used for reducing the bend. Equivalent to the diameter of a circle to remove

     Returns
     -------
     list of QgsFeatures
        List of reduces QgsFeatures
     """

    rb_return = RbResults()
    
    # Create the list of RbLineString and RbPoint to process
    rb_features = _create_rb_feature(qgs_in_features)
    rb_return.in_nbr_features = len(qgs_in_features)
    print ("test:", rb_return.in_nbr_features)

    # Create the list of RbGeom ==> List of geometry to reduce the bend
    rb_geoms = []
    for rb_feature in rb_features:
        rb_geoms += rb_feature.get_rb_geom()

    # Create the RbCollection a spatial index to accelerate search
    rb_collection = RbCollection()
    rb_collection.add_features(rb_geoms)

    # Execute the bend reduction
    _manage_reduce_bend(rb_geoms, rb_collection, diameter_tol)

    # Recreate the QgsFeature
    qgs_features_out = [rb_feature.get_qgs_feature() for rb_feature in rb_features]

    # Calculate return values
    print (rb_return.in_nbr_features)
    rb_return.out_nbr_features = len(qgs_features_out)
    rb_return.nbr_bends_reduced = sum([rb_geom.nbr_bend_reduced for rb_geom in rb_geoms])
    rb_return.nbr_bends_detected = sum([rb_geom.nbr_bend_detected for rb_geom in rb_geoms])
    rb_return.qgs_features_out = qgs_features_out

    return rb_return
