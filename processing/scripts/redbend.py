import os, sys
file = os.path.abspath(__file__)
pathname = os.path.dirname(os.path.abspath(__file__))
sys.path.append(pathname)
print ("remove comment...")
###from algo_reduce_bend import reduce_bends, is_line_string, is_polygon

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













# Remaining modifications:
#  - attribute for isClosed instead off calling each time isClosed()
#  - put comment in the code
#  - Manage case where a closed line string is of the same orientation except the start/end which is different
#  - Have a dynamic epsilon for some cases where it's needed
#  - Add some comprehension list instead of for loop
#  - edit line with a smooth line instead of a stragth line
#  - reorient the bend solely on the second pass as it will mobe start/end probably on a greater bend
#  - in the bend flagging process prioritize the bend that goes outside the polygon first (for polygon)
#  - Test reduce_bend performance with the profiler

from abc import ABC, abstractmethod
import math
from qgis.core import QgsLineString, QgsPolygon, QgsWkbTypes, QgsSpatialIndex, \
    QgsGeometry, QgsGeometryUtils, QgsVertexId

# Define global constant
GeoSim_EPSILON = 1.0E-6
GeoSim_CW = 0
GeoSim_ACW = -1


def is_polygon(feature_type):
    if feature_type in [QgsWkbTypes.Polygon, QgsWkbTypes.Polygon25D, QgsWkbTypes.PolygonZ, QgsWkbTypes.PolygonM,
                        QgsWkbTypes.PolygonZM]:
        val = True
    else:
        val = False

    return val


def is_line_string(feature_type):
    """Routine to test if the feature type is a line string

    Parameters
    ----------
    feature_type : int
        Feature type to test

    Returns
    -------
    Boolean
        True: it's a line string; False: Otherwise
    """

    if feature_type in [QgsWkbTypes.LineString, QgsWkbTypes.LineString25D, QgsWkbTypes.LineStringZ,
                        QgsWkbTypes.LineStringM, QgsWkbTypes.LineStringZM]:
        val = True
    else:
        val = False

    return val


def is_point(feature_type):
    if feature_type in [QgsWkbTypes.Point, QgsWkbTypes.Point25D, QgsWkbTypes.PointM, QgsWkbTypes.PointZ,
                        QgsWkbTypes.PointZM]:
        val = True
    else:
        val = False

    return val


def reduce_bends(qgs_in_features, diameter_tol, feedback=None, flag_del_outer=False, flag_del_inner=False):
    """Main routine for bend reduction

    Parameters
    ----------
    qgs_in_features : List of QgsFeatures
         List of QgsFeatures to reduce
    diameter_tol : Float
         Tolerance used for reducing the bend. Equivalent to the diameter of a circle to remove
    flag_del_outer : Boolean (default is False)
         Flag indicating if we delete feature that have an exterior ring area below the tolerance
    flag_del_inner : Boolean (default is False)
         Flag indicating if we delete inner rings that have an area below the tolerance

     Returns
     -------
     list of QgsFeatures
        List of reduces QgsFeatures
     """

    import cProfile, pstats, io
    from pstats import SortKey
    pr = cProfile.Profile()
    pr.enable()

    rb_results = RbResults()

    # Create the list of RbLineString and RbPoint to process
    rb_features = _create_rb_feature(qgs_in_features)
    rb_results.in_nbr_features = len(qgs_in_features)

    # Delete the outer or inner ring below the diameter tolerance
    if flag_del_outer or flag_del_inner:
        del_outer_inner_ring(rb_features, rb_results, diameter_tol, flag_del_outer, flag_del_inner)

    # Create the list of RbGeom ==> List of geometry to reduce the bend
    rb_geoms = []
    for rb_feature in rb_features:
        rb_geoms += rb_feature.get_rb_geom()

    # Create the RbCollection a spatial index to accelerate search
    rb_collection = RbCollection()
    rb_collection.add_features(rb_geoms)

    # Execute the bend reduction
    _manage_reduce_bend(rb_geoms, rb_collection, rb_results, diameter_tol, feedback)

    # Recreate the QgsFeature
    qgs_features_out = [rb_feature.get_qgs_feature() for rb_feature in rb_features]

    # Calculate return values
    rb_results.out_nbr_features = len(qgs_features_out)
    rb_results.nbr_bends_reduced = sum([rb_geom.nbr_bend_reduced for rb_geom in rb_geoms])
    rb_results.nbr_bends_detected = sum([rb_geom.nbr_bend_detected for rb_geom in rb_geoms])
    rb_results.qgs_features_out = qgs_features_out

    pr.disable()
    s = io.StringIO()
    sortby = SortKey.CUMULATIVE
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print(s.getvalue())

    return rb_results


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
        self.rb_geom[0].is_simplest = True  # A point cannot be reduced
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
        RbGeom._id_counter += 1
        self.qgs_geom = QgsGeometry(qgs_geom)
        self.qgs_rectangle = self.qgs_geom.boundingBox()
        self.bends = []
        RbGeom._id_counter += 1
        self.nbr_bend_reduced = 0
        if is_point(qgs_geom.wkbType()):
            # A point cannot be reduced
            self.is_simplest = True
            self.nbr_bend_detected = 0
        else:
            self.is_simplest = False
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
                    new_start_end = (bend_area_ok.j + bend_area_ok.i) // 2
                    new_qgs_points = qgs_points[new_start_end:] + qgs_points[1:new_start_end + 1]
                    self.qgs_geom = QgsGeometry(QgsLineString(new_qgs_points))
                    detect_bends(self)


class Bend:

    def __init__(self, i, j, qgs_points):

        self.i = i
        self.j = j
        qgs_pol_bend = QgsPolygon(QgsLineString(qgs_points[i:j + 1]))
        self.bend_area = qgs_pol_bend.area()
        self.bend_perimeter = qgs_pol_bend.perimeter()
        self.qgs_geom_bend = QgsGeometry(qgs_pol_bend.clone())
        self.adj_area = calculate_adj_area(self.bend_area, self.bend_perimeter)
        self.to_reduce = False
        self.qgs_geom_new_subline = None
        self.qgs_geom_new_subline_trimmed = None

    def reduce(self, rb_geom):

        for ind in range(self.j - 1, self.i, -1):  # Process the vertice in reverse order
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
        self.in_nbr_features = None
        self.out_nbr_features = None
        self.nbr_bend_reduced = None
        self.qgs_features_out = None
        self.nbr_hole_del = 0
        self.nbr_pol_del = 0

    # def dummy():
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
        qgs_geom = qgs_feature.geometry()  # extract the Geometry
        feature_type = qgs_geom.wkbType()

        if is_polygon(feature_type):
            rb_features.append(RbPolygon(qgs_feature))
        elif is_line_string(feature_type):
            rb_features.append(RbLineString(qgs_feature))
        elif is_point(feature_type):
            rb_features.append(RbPoint(qgs_feature))
        else:
            raise Exception("Unsupported GeometryType: {}".format(qgs_geom.wkbType()))

    return rb_features


def angle_between_three_points(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y):

    angle1 = math.atan2(p0_y - p1_y, p0_x - p1_x )
    angle2 = math.atan2(p2_y - p1_y, p2_x - p1_x)

    angle = angle1 - angle2

    # Normalizing angle
    clippedAngle = angle
    if clippedAngle >= math.pi * 2 or clippedAngle <= -2 * math.pi:
        clippedAngle = clippedAngle % 2*math.pi
    if clippedAngle < 0.0:
        clippedAngle += 2 * math.pi

    return clippedAngle

#def create_coords_xy(qgs_geom):

#    str_wkt = qgs_geom.asWkt()  # Convert to WKT
#    start = str_wkt.find("(")  # Find the position of the first parenthesis
#    end = str_wkt.find(")")  # Find the position of the first parenthesis
#    str_wkt = (str_wkt[start + 1:end])  # Trim the line string
#    lst_coord_str = str_wkt.split(", ")
#    lst_coord_xy = []
#    for coord_str in lst_coord_str:
#        coord_xy = (coord_str.split(" "))
#        x = float(coord_xy[0])
#        y = float(coord_xy[1])
#        lst_coord_xy.append((x, y))
#
#    return lst_coord_xy


def pseudo_reduce(qgs_line_string):
    # Delete almost duplicate point
#    qgs_line_string = rb_geom.qgs_geom.constGet()
    qgs_line_string.removeDuplicateNodes(GeoSim_EPSILON)

    # Delete co-linear points on the line angle of near 0 or 180 degrees
#    qgs_points = qgs_line_string.points()
#    coords_xy = create_coords_xy(rb_geom.qgs_geom)
    coords_xy = [(qgs_point.x(),qgs_point.y()) for qgs_point in qgs_line_string.points()]
    num_points = len(coords_xy)
#    num_points = len(coords_xy)
    i = 1
    while i <= num_points - 2 and num_points >= 3:
#        p0_x, p0_y = qgs_points[i - 1].x(), qgs_points[i - 1].y()
#        p1_x, p1_y = qgs_points[i].x(), qgs_points[i].y()
#        p2_x, p2_y = qgs_points[i + 1].x(), qgs_points[i + 1].y()
        p0_x, p0_y = coords_xy[i - 1][0], coords_xy[i - 1][1]
        p1_x, p1_y = coords_xy[i][0], coords_xy[i][1]
        p2_x, p2_y = coords_xy[i + 1][0], coords_xy[i + 1][1]

        angle = QgsGeometryUtils.angleBetweenThreePoints(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y)
#        angle = angle_between_three_points(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y)
        #        if abs(angle-angle_prime) > GeoSim_EPSILON:
        #            print ("calcul angle pas bon...")
        #            0/0
        if abs(angle - math.pi) <= GeoSim_EPSILON or abs(angle) <= GeoSim_EPSILON:
            # Co-linear point or flat angle delete the current point
            #del coords_xy[i]
            del coords_xy[i]
            vertex_id = QgsVertexId(_part=0, _ring=0, _vertex=i, _type=QgsVertexId.SegmentVertex)
            qgs_line_string.deleteVertex(vertex_id)
            num_points -= 1
        else:
            i += 1

    return coords_xy


def detect_bends(rb_geom):

    qgs_line_string = rb_geom.qgs_geom.constGet()

    # Delete co-linear and almost duplicate point
    coords_xy = pseudo_reduce(qgs_line_string)

    qgs_points = qgs_line_string.points()
##    coords_xy = [(qgs_point.x(), qgs_point.y()) for qgs_point in qgs_points]

    if len(qgs_points) != len(coords_xy):
        raise Exception ("Internal corruption detected in module: detect_bends")

#    coords_xy = create_coords_xy(rb_geom.qgs_geom)

    num_points = len(coords_xy)
    i = 1
    angles = []
    # Extract the angles at each vertice except start and end vertice
    while i <= num_points - 2 and num_points >= 3:
        p0_x, p0_y = coords_xy[i - 1][0], coords_xy[i - 1][1]
        p1_x, p1_y = coords_xy[i][0], coords_xy[i][1]
        p2_x, p2_y = coords_xy[i + 1][0], coords_xy[i + 1][1]
        angles.append(QgsGeometryUtils.angleBetweenThreePoints(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y))
        #        angle_prime = angle_between_three_points(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y)
        #        if abs(angles[-1] - angle_prime) > GeoSim_EPSILON:
        #            print("calcul angle pas bon...")
        #            0 / 0
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
                pass  # Nothing to do
            else:
                # Change of direction.  Create a bend
                rb_geom.bends.append(Bend(start, i + 1, qgs_points))
                direction = orientation[i]
                start = i
            i += 1  # Go look for the next vertice

        rb_geom.bends.append(Bend(start, i + 1, qgs_points))  # Create the last bend

        if qgs_line_string.isClosed() and \
                i + 2 - start == num_points:
            # A closed lines with only one bend cannot be reduced
            rb_geom.bends = []  # Reset the bends

    if len(rb_geom.bends) == 0:
        # There is no bend so nothing to simplify
        rb_geom.is_simplest = True

    if rb_geom.nbr_bend_detected is None:
        # For the first pass set the number of bend detected (for statistics purpose)
        rb_geom.nbr_bend_detected = len(rb_geom.bends)

    return


#def is_bend_reduction_terminated(last_nbr_bend_reduced, rb_geoms):
#    new_nbr_bend_reduced = sum(rb_geom.nbr_bend_reduced for rb_geom in rb_geoms)
#    if new_nbr_bend_reduced == last_nbr_bend_reduced:
#        bend_to_reduce = False
#    else:
#        bend_to_reduce = True
#
#    return bend_to_reduce


def remove_rb_geoms_done(rb_geoms, rb_geoms_done):
    for i in reversed(range(len(rb_geoms))):
        if rb_geoms[i].is_simplest:
            rb_geoms_done.append(rb_geoms[i])
            del rb_geoms[i]

    #    for i, rb_geom in enumerate(rb_geoms):
    #        if rb_geom.is_simplest:
    #            to_move.append(i)

    #    # Move from rb_geoms to rb_geoms_done
    #    print ("to_move len: ", len(to_move))
    #    for i in reversed(range(len(to_move))):   # Process the list in reverse order
    #        rb_geoms_done.append(rb_geoms[i])
    #        del rb_geoms[i]

    return


def calculate_adj_area(area, perimeter):
    compactness_index = 4 * area * math.pi / perimeter ** 2
    adj_area = area * (.75 / compactness_index)

    return adj_area


def calculate_min_adj_area(diameter_tol):
    min_adj_area = .75 * math.pi * (diameter_tol / 2.) ** 2

    return min_adj_area


def flag_bend_to_reduce(rb_geom, diameter_tol):
    # Minimum adjusted area used to find bend to reduce
    min_adj_area = calculate_min_adj_area(diameter_tol)

    if rb_geom.qgs_geom.constGet().isClosed and len(rb_geom.bends) >= 3:
        # The closed line has been rotated and the start/end point lie on a bend that do not need to be reduced
        del rb_geom.bends[0]  # Remove the first bend
        del rb_geom.bends[-1]  # Remove the last bend

    #    lst_bends = []
    #    for i, bend in enumerate(rb_geom.bends):
    #        lst_bends.append((bend.adj_area, i))
    lst_bends = [(bend.adj_area, i) for i, bend in enumerate(rb_geom.bends)]

    # Sort from smallest to biggest bend area
    lst_bends.sort(key=lambda item: item[0])
    #    print ("min_adj_area", min_adj_area)
    #    print (lst_bends)
    start = 0
    end = len(rb_geom.bends) - 1

    for (adj_area, i) in lst_bends:
        if adj_area <= min_adj_area:
            if len(lst_bends) == 1:
                rb_geom.bends[i].to_reduce = True  # Only one bend process it...
            else:
                if i == start:
                    if rb_geom.bends[i + 1].to_reduce:
                        pass  # Cannot reduce two bend adjacent
                    else:
                        rb_geom.bends[i].to_reduce = True
                elif i == end:
                    if rb_geom.bends[i - 1].to_reduce:
                        pass  # Cannot reduce two bend adjacent
                    else:
                        rb_geom.bends[i].to_reduce = True
                elif rb_geom.bends[i - 1].to_reduce or rb_geom.bends[i + 1].to_reduce:
                    pass  # Cannot reduce two bend adjacent
                else:
                    rb_geom.bends[i].to_reduce = True
        else:
            # Over minimum adjusted area
            break

    if len(lst_bends) == 0 or lst_bends[0][0] >= min_adj_area:
        rb_geom.is_simplest = True


def validate_spatial_constraints(validate_is_simple, ind, rb_geom, rb_collection):

    check_constraints = True
    bend = rb_geom.bends[ind]
    qgs_geom_line_string = rb_geom.qgs_geom 

    # First: check if the bend reduce line string is an OGC simple line
    # We test with a tiny smaller line to ease the testing and false positive error
    if validate_is_simple:
        qgs_geom_new_subline_trimmed = bend.get_new_subline_trimmed(rb_geom)
        qgs_geom_new_sub_trim_engine = QgsGeometry.createGeometryEngine(qgs_geom_new_subline_trimmed.constGet())
        if qgs_geom_new_sub_trim_engine.disjoint(qgs_geom_line_string.constGet()):
            # Everything is OK
            pass
        else:
            # The new sub line intersect the line itself. The result would create a non OGC simple line
            check_constraints = False
            print (qgs_geom_new_sub_trim_engine.intersection(qgs_geom_line_string.constGet()))

    # Second: check that the new line does not intersect any other line or points
    if check_constraints:
        qgs_rectangle = bend.qgs_geom_bend.boundingBox()
        qgs_geom_potentials = rb_collection.get_features(qgs_rectangle, [rb_geom.id])
        qgs_geom_new_subline = bend.get_new_subline(rb_geom)
        qgs_geom_engine_new_subline = QgsGeometry.createGeometryEngine(qgs_geom_new_subline.constGet())
        for qgs_geom_potential in qgs_geom_potentials:
#            if not qgs_geom_new_subline.disjoint(qgs_geom_potential):
            if qgs_geom_engine_new_subline.intersects(qgs_geom_potential.constGet()):
                # The bend area intersects with a point
                check_constraints = False
                break

    # Third: check that inside the bend to reduce there is no feature completely inside it.  This would cause a
    # sidedness or relative position error
    if check_constraints:
#        qgs_geom_engine_bend_area = QgsGeometry.createGeometryEngine(bend.qgs_geom_bend)
#        qgs_geom_bend = QgsGeometry(bend.qgs_geom_bend.clone())
        for qgs_geom_potential in qgs_geom_potentials:
#            if qgs_geom_bend_area.contains(qgs_geom_potential.constGet()):
            if bend.qgs_geom_bend.contains(qgs_geom_potential):
                # A feature is totally located inside
                check_constraints = False
                break

    return check_constraints


def process_bends(rb_geom, rb_collection):

    for pass_num in (0,1):
        nbr_bend_reduced = 0
        if pass_num == 0:
            if rb_geom.qgs_geom.isSimple():
                qgs_points_copy = rb_geom.qgs_geom.constGet().points()
                validate_is_simple = False
            else:
                # The line is already not simple.  Should not happen but algothim managed this case
                validate_is_simple = True
                print("OGC non valide avant pass 0")
        else:
            # A second pass indicates a problem of line self intersecting (manage intersection of each bend reduction)
            validate_is_simple = True

        for ind in reversed(range(len(rb_geom.bends))):
            bend = rb_geom.bends[ind]
            if bend.to_reduce:
                # Check spatial constraints
                spatial_constraints = validate_spatial_constraints(validate_is_simple, ind, rb_geom, rb_collection)
                if spatial_constraints:
                    bend.reduce(rb_geom)
                    nbr_bend_reduced += 1

        if pass_num == 0:
            # Validate for line self intersection
            if rb_geom.qgs_geom.isSimple():
                # The second pass is not necessary
                break
            else:
                # Undo the edit and recreate the original QgsLineString
                qgs_line_string = QgsLineString(qgs_points_copy)
                rb_geom.qgs_geom = QgsGeometry(qgs_line_string.clone())
                print("OGC non valide apres pass 0")
        else:
            # Nothing to validate in a second pass
            pass

    return nbr_bend_reduced



def _manage_reduce_bend(rb_geoms, rb_collection, rb_results, diameter_tol, feedback):

    rb_geoms_done = []
    nbr_pass = 0
    previous_pass_nbr_bends = -1
    current_pass_nbr_bends = 0
    nbr_geoms = 100.0 / len(rb_geoms) if len(rb_geoms) >= 1 else 0
    while previous_pass_nbr_bends != current_pass_nbr_bends:
        remove_rb_geoms_done(rb_geoms, rb_geoms_done)  # Remove feature done to accelerate process
        # set the progress bar
        if feedback is not None:
            if len(rb_geoms_done) == 0:
                feedback.setProgress(1)
            else:
                feedback.setProgress(int(len(rb_geoms_done) * nbr_geoms))
        previous_pass_nbr_bends = current_pass_nbr_bends
        current_pass_nbr_bends = 0
        for rb_geom in rb_geoms:
            detect_bends(rb_geom)
            if nbr_pass == 0:
                # Edit the start/end point for closed QgsLineString
                rb_geom.edit_closed_line(diameter_tol)
            flag_bend_to_reduce(rb_geom, diameter_tol)
            current_pass_nbr_bends += process_bends(rb_geom, rb_collection)
        # Check if all bend are processed
#        is_terminated = is_bend_reduction_terminated(last_nbr_bend_reduced, rb_geoms)
        nbr_pass += 1

    # Reset the rb_geoms list
    rb_geoms += rb_geoms_done
    rb_results.nbr_pass = nbr_pass

    return


def extract_polygon_attributes(qgs_geom):
    qgs_line_string = qgs_geom.constGet()
    qgs_pol = QgsPolygon(qgs_line_string.clone())
    area = qgs_pol.area()
    perimeter = qgs_pol.perimeter()

    return (area, perimeter)


def del_outer_inner_ring(rb_features, rb_results, diameter_tol, flag_del_outer, flag_del_inner):
    # Loop over each rb_features
    for i in reversed(range(len(rb_features))):
        if isinstance(rb_features[i], RbPolygon):  # Only process Polygon
            min_adj_area = calculate_min_adj_area(diameter_tol)
            for j in reversed(range(len(rb_features[i].rb_geom))):
                area, perimeter = extract_polygon_attributes(rb_features[i].rb_geom[j].qgs_geom)
                adj_area = calculate_adj_area(area, perimeter)
                if j == 0:
                    # Process the exterior ring
                    if flag_del_outer:
                        if adj_area < min_adj_area:
                            del rb_features[i]  # Delete the rb_feature
                            rb_results.nbr_pol_del += 1
                            break
                else:
                    # Process an interior ring
                    if flag_del_inner:
                        if adj_area < min_adj_area:
                            del rb_features[i].rb_geom[j]  # Delete the ring
                            rb_results.nbr_hole_del += 1
