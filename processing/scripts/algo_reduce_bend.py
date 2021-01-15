#Remaining modifications:
#  - put comment in the code
#  - verify if progrssive erosion bring a plus (min_pass)
#  - Add some comprehension list instead of for loop
#  - edit line with a smooth line instead of a straight line
#  - in the bend flagging process prioritize the bend that goes outside the polygon first (for polygon)
#  - put some code to correct J-Bend

from abc import ABC, abstractmethod
import sys, math
from qgis.core import QgsFeature, QgsPoint, QgsPointXY, QgsLineString, QgsPolygon, QgsWkbTypes, QgsSpatialIndex, \
    QgsGeometry, QgsGeometryUtils, QgsRectangle, QgsProcessingException

# Define global constant
ANTI_CLOCK_WISE = -1
CLOCK_WISE = 0


class RbFeature(ABC):
    _id_counter = 0

    @staticmethod
    def is_point(feature_type):
        if feature_type in [QgsWkbTypes.Point, QgsWkbTypes.Point25D, QgsWkbTypes.PointM, QgsWkbTypes.PointZ,
                            QgsWkbTypes.PointZM]:
            val = True
        else:
            val = False

        return val

    @staticmethod
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

    @staticmethod
    def is_polygon(feature_type):
        if feature_type in [QgsWkbTypes.Polygon, QgsWkbTypes.Polygon25D, QgsWkbTypes.PolygonZ, QgsWkbTypes.PolygonM,
                            QgsWkbTypes.PolygonZM]:
            val = True
        else:
            val = False

        return val

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
        self.rb_geom = [RbGeom(qgs_geom, QgsWkbTypes.Polygon) for qgs_geom in self.qgs_geom]
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

        self.rb_geom = [RbGeom(self.qgs_geom, QgsWkbTypes.LineString)]
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

        self.rb_geom = [RbGeom(self.qgs_geom, QgsWkbTypes.Point)]
        self.rb_geom[0].is_simplest = True  # A point cannot be reduced
        self.qgs_geom = None

    def get_rb_geom(self):
        return self.rb_geom

    def get_qgs_feature(self):
        qgs_geom = QgsGeometry(self.rb_geom[0].qgs_geom.constGet().clone())
        self.qgs_feature.setGeometry(qgs_geom)
        return self.qgs_feature


class GeomLineSegment:

    def __init__(self):

        self._dict_index = {}  # Dictionary of spatial index for line segment per
        self._id = 0

    def _get_next_id(self):

        self._id += 1

        return self._id

    def add_geom(self, rb_geom):

        if rb_geom.qgs_geom.wkbType() == QgsWkbTypes.LineString:
            spatial_index = QgsSpatialIndex(flags=QgsSpatialIndex.FlagStoreFeatureGeometries)
            qgs_points = rb_geom.qgs_geom.constGet().points()
            num_points = len(qgs_points)
            qgs_features = []
#            if rb_geom.is_closed:
#                nbr_segment = num_points  # Take into account that the start/end point (circular array)
#            else:
            nbr_segment = num_points-1
            for i in range(nbr_segment):
                qgs_feature= QgsFeature(id=self._get_next_id())
                qgs_feature.setGeometry(QgsLineString(qgs_points[i], qgs_points[(i+1)%num_points]))
                qgs_features.append(qgs_feature)

            # Add the line segment in the spatial index
            spatial_index.addFeatures(qgs_features)
            self._dict_index[rb_geom.id] = spatial_index

        return

    def delete_line_segment(self, qgs_geom_id, qgs_pnt0, qgs_pnt1):

        spatial_index = self._dict_index[qgs_geom_id]
        qgs_mid_point = QgsGeometryUtils.midpoint(qgs_pnt0, qgs_pnt1)
        qgs_rectangle = qgs_mid_point.boundingBox()
        qgs_rectangle.grow(ReduceBend.ZERO_RELATIVE*100)
        ids = spatial_index.intersects(qgs_rectangle)
        for id in ids:
            qgs_geom_line = spatial_index.geometry(id)  # Extract geometry from index
            qgs_pnt_start = qgs_geom_line.vertexAt(0)
            qgs_pnt_end = qgs_geom_line.vertexAt(1)
            if qgs_pnt_start.distance(qgs_pnt0) <= ReduceBend.ZERO_RELATIVE and \
               qgs_pnt_end.distance(qgs_pnt1) <= ReduceBend.ZERO_RELATIVE:
                feature = QgsFeature(id=id)
                feature.setGeometry(QgsLineString([qgs_pnt_start, qgs_pnt_end]))
                if (spatial_index.deleteFeature(feature)):
                    deleted = True
                    break
                else:
                    raise Exception (QgsProcessingException("Unable to delete entry in QgsSpatialIndex..."))
            else:
                deleted = False

        if not deleted:
            raise Exception(QgsProcessingException("Internal structure corruption..."))

        return

    def add_line_segment(self, rb_geom_id, qgs_pnt0, qgs_pnt1):

        feature = QgsFeature(id=self._get_next_id())
        feature.setGeometry(QgsLineString(qgs_pnt0, qgs_pnt1))
        self._dict_index[rb_geom_id].addFeature(feature)

        return

    def get_line_segments(self, geom_id, qgs_rectangle):

        spatial_index = self._dict_index[geom_id]  # Extract the appropriate spatial index
        ids = spatial_index.intersects(qgs_rectangle)
        # Transform the keys into rb_features
        qgs_geom_line_segments = [spatial_index.geometry(id) for id in ids]

        return qgs_geom_line_segments


class RbCollection(object):

    __slots__ = ('_spatial_index', 'space_size', '_dict_qgs_geom', '_dict_geom_b_boxes', '_id_b_box',
                 '_geom_line_segment', 'rb_results')

    def __init__(self, rb_results, nbr_sub_space):
        self._spatial_index = QgsSpatialIndex(flags=QgsSpatialIndex.FlagStoreFeatureGeometries)
        self.space_size = ReduceBend.MAP_RANGE / nbr_sub_space
        self._dict_qgs_geom = {}
        self._dict_geom_b_boxes = {}
        self._id_b_box = 0
        self._geom_line_segment = GeomLineSegment()
        self.rb_results = rb_results

        return

    def _get_id_b_box(self):

        self._id_b_box += 1

        return self._id_b_box

    def _get_b_box_xy(self, qgs_b_box):

        x_min = qgs_b_box.xMinimum()
        y_min = qgs_b_box.yMinimum()
        x_max = qgs_b_box.xMaximum()
        y_max = qgs_b_box.yMaximum()

        return x_min, y_min, x_max, y_max

    def _polygon_from_rectangle(self, qgs_b_box):

        x_min, y_min, x_max, y_max = self._get_b_box_xy(qgs_b_box)
        qgs_pol = QgsPolygon(QgsLineString((QgsPoint(x_min, y_min), QgsPoint(x_max, y_min),
                                            QgsPoint(x_max, y_max), QgsPoint(x_min, y_max))))

        return qgs_pol

    def _subdivide_b_box(self, qgs_geom, qgs_b_box, lst_b_box):

        if max(qgs_b_box.width(),qgs_b_box.height()) <= self.space_size:
            # No bounding box should be 0 or near 0 in width or in height
            # Also avoid some side effects when bends are parallel to the boundind box
            qgs_b_box.grow(ReduceBend.ZERO_RELATIVE*100)
            if qgs_geom.intersects(qgs_b_box):
#                x_min, y_min, x_max, y_max = self._get_b_box_xy(qgs_b_box)
                qgs_feature_b_box = QgsFeature(id=self._get_id_b_box())
                qgs_geom_b_box = self._polygon_from_rectangle(qgs_b_box)
                qgs_feature_b_box.setGeometry(qgs_geom_b_box)
                lst_b_box += [qgs_feature_b_box]  # Keep that bounding box
        else:
            # Sub divide this bounding box
            x_min, y_min, x_max, y_max = self._get_b_box_xy(qgs_b_box)
            x_mid = (x_min + x_max) / 2.
            y_mid = (y_min + y_max) / 2.
            if x_max-x_min < self.space_size:
                # Only create north-south quadrant
                qgs_rect_s = QgsRectangle(x_min, y_min, x_max, y_mid)
                qgs_rect_n = QgsRectangle(x_min, y_mid, x_max, y_max)
                qgs_rectangles = [qgs_rect_s, qgs_rect_n]
            elif y_max-y_min < self.space_size:
                # Only create east-west quadrant
                qgs_rect_w = QgsRectangle(x_min, y_min, x_mid, y_max)
                qgs_rect_e = QgsRectangle(x_mid, y_min, x_max, y_max)
                qgs_rectangles = [qgs_rect_w, qgs_rect_e]
            else:
                # Split the current bounding box into four smaller quadrant (bounding boxes)
                qgs_rect_sw = QgsRectangle(x_min, y_min, x_mid, y_mid)
                qgs_rect_se = QgsRectangle(x_mid, y_min, x_max, y_mid)
                qgs_rect_nw = QgsRectangle(x_min, y_mid, x_mid, y_max)
                qgs_rect_ne = QgsRectangle(x_mid, y_mid, x_max, y_max)
                qgs_rectangles = [qgs_rect_sw, qgs_rect_se, qgs_rect_nw, qgs_rect_ne]
            # Loop over each quadrant
            for qgs_rect in qgs_rectangles:
                lst_b_box = self._subdivide_b_box(qgs_geom, qgs_rect, lst_b_box)

        return lst_b_box


    def create_index_b_boxes(self, rb_geom_id, qgs_geom):

        qgs_features_b_box = self._subdivide_b_box(qgs_geom, qgs_geom.boundingBox(), [])
        b_box_ids = [qgs_feature_b_box.id() for qgs_feature_b_box in qgs_features_b_box]
        for b_box_id in b_box_ids:
            self._dict_geom_b_boxes[b_box_id] = rb_geom_id  # Keep link between b_box id and the rb_geom id

        self._spatial_index.addFeatures(qgs_features_b_box)  # Add the bounding boxes into the index

        return len(b_box_ids)


    def add_features(self, rb_geoms):

        b_boxes = 0
        for rb_geom in rb_geoms:
            b_boxes += self.create_index_b_boxes(rb_geom.id, rb_geom.qgs_geom)
            self._dict_qgs_geom[rb_geom.id] = rb_geom.qgs_geom  # Keep a reference to the qgs_geom
            self._geom_line_segment.add_geom(rb_geom)  #  Create the line segment index for the geometry

        self.rb_results.lines_log_info.append("Number of bounding box created: {}".format(b_boxes))

        return


    def get_features(self, qgs_rectangle, drop_ids=[]):

        b_box_ids = self._spatial_index.intersects(qgs_rectangle)
        # From the b_box_ids get the qgs_geom ids
        geom_ids = [self._dict_geom_b_boxes[b_box_id] for b_box_id in b_box_ids]
        geom_ids = list(set(geom_ids))  # Remove duplicate geom_ids in the list
        # From the id extract the qgs_geom
        qgs_geom = [self._dict_qgs_geom[geom_id] for geom_id in geom_ids if geom_id not in drop_ids]

        return qgs_geom


    def delete_vertex(self, rb_geom, v_ids_to_del):

        qgs_geom = rb_geom.qgs_geom
        v_ids_to_process = [v_ids_to_del[0]-1] + v_ids_to_del + [v_ids_to_del[-1]+1]
        for i in range(len(v_ids_to_process)-1):
            qgs_pnt_first = qgs_geom.vertexAt(v_ids_to_process[i])
            qgs_pnt_last = qgs_geom.vertexAt(v_ids_to_process[i+1])
            self._geom_line_segment.delete_line_segment(rb_geom.id, qgs_pnt_first, qgs_pnt_last)

        # Add the new segment line
        qgs_pnt_0 = qgs_geom.vertexAt(v_ids_to_process[0])
        qgs_pnt_1 = qgs_geom.vertexAt(v_ids_to_process[-1])
        self._geom_line_segment.add_line_segment(rb_geom.id, qgs_pnt_0, qgs_pnt_1)

        # Check that the new line segment added is completely covered by the bounding boxes of the geometry
        qgs_rectangle = QgsRectangle(QgsPointXY(qgs_pnt_0), QgsPointXY(qgs_pnt_1))
        qgs_geom_line_string = QgsGeometry(QgsLineString(qgs_pnt_0, qgs_pnt_1))
        b_box_ids = self._spatial_index.intersects(qgs_rectangle)
        for b_box_id in b_box_ids:
            if self._dict_geom_b_boxes[b_box_id] == rb_geom.id:
                qgs_pol = self._spatial_index.geometry(id=b_box_id)
                if qgs_geom_line_string.within(qgs_pol):
                    line_segment_covered = True
                    break
                else:
                    line_segment_covered = False
        if (not line_segment_covered):
            self.create_index_b_boxes(rb_geom.id, qgs_geom_line_string)

        # Now that the line segment spatial structure is updated let's delete the vertex
        for v_id_to_del in reversed(v_ids_to_del):
            qgs_geom.deleteVertex(v_id_to_del)

        return


    def get_line_segments(self, geom_id, qgs_rectangle):

        qgs_geom_line_segments = self._geom_line_segment.get_line_segments(geom_id, qgs_rectangle)

        return qgs_geom_line_segments


class RbGeom:

    __slots__ = ('id', 'original_geom_type', 'is_simplest', 'qgs_geom', 'qgs_rectangle', 'bends', 'nbr_bend_reduced',
                 'is_closed', 'need_pivot')

    _id_counter = 0

    def __init__(self, qgs_abs_geom, original_geom_type):

        self.id = self.next_id()
        self.original_geom_type = original_geom_type
        qgs_geometry = qgs_abs_geom.constGet()
        self.qgs_geom = QgsGeometry(qgs_geometry.clone())
        self.is_simplest = False
        self.is_closed = False
        self.need_pivot = False
        if self.original_geom_type == QgsWkbTypes.Point:
            self.is_simplest = True  # A point cannot be simplified
        elif self.original_geom_type == QgsWkbTypes.LineString:
            if qgs_geometry.length() <= ReduceBend.ZERO_RELATIVE:
                self.is_simplest = True  # Zero length line. Do not try to simplify
            else:
                if qgs_geometry.isClosed():
                    if abs(qgs_geometry.sumUpArea()) > ReduceBend.ZERO_RELATIVE:
                        self.is_closed = True
                        self.need_pivot = True
#                        self.orientation = qgs_geometry.orientation()
                    else:
                        self.is_simplest = True # Zero area polygon,  Do not try to simplify
        else:
            # It's a polygon
            if abs(qgs_geometry.sumUpArea()) > ReduceBend.ZERO_RELATIVE:
                self.is_closed = True
                self.need_pivot = True
            else:
                self.is_simplest = True  # Zero area polygon. Do not simplify the closed line

        self.qgs_rectangle = self.qgs_geom.boundingBox()
        self.bends = []
        self.nbr_bend_reduced = 0

    def next_id(self):

        RbGeom._id_counter += 1

        return (RbGeom._id_counter)

    def get_angles(self):

        if self.original_geom_type != QgsWkbTypes.Point:
            qgs_line_string = self.qgs_geom.constGet()
            num_xy = qgs_line_string.numPoints()
            xy = [(qgs_line_string.xAt(i), qgs_line_string.yAt(i)) for i in range(num_xy)]
            if self.is_closed:
                # Add two vertice at the start/end for the circularity of a closed line
                end_xy = xy[-2]  # Not that last because it's the same position as the first vertice
                xy.insert(0, end_xy)
            # Extract the angles at each vertice except start and end vertice
            angles = [QgsGeometryUtils.angleBetweenThreePoints(xy[i-1][0], xy[i-1][1], xy[i][0], xy[i][1],
                                                               xy[i+1][0], xy[i+1][1]) for i in range(1, len(xy)-1)]

        return angles


class Bend:

    def __init__(self, i, j, qgs_polygon):

        __slots__ = ('i', 'j', 'area', 'perimeter', 'adj_area', 'to_reduce', 'qgs_geom_new_subline',
                     'qgs_geom_new_subline_trimmed', 'orientation')

        self.i = i
        self.j = j
        self.qgs_geom_bend = QgsGeometry(qgs_polygon.clone())
        self.area = self.qgs_geom_bend.area()
        self.perimeter = self.qgs_geom_bend.length()
        self.adj_area = ReduceBend.calculate_adj_area(self.area, self.perimeter)
        self.to_reduce = False
        self.orientation = None
        self.qgs_geom_new_subline = None
        self.qgs_geom_new_subline_trimmed = None

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
            if qgs_ls_new_line.length() >= ReduceBend.ZERO_RELATIVE*100.:
                qgs_pnt_i_trimmed = qgs_ls_new_line.interpolatePoint(ReduceBend.ZERO_RELATIVE)
                qgs_pnt_j_trimmed = qgs_ls_new_line.interpolatePoint(qgs_ls_new_line.length() - ReduceBend.ZERO_RELATIVE)
                qgs_ls_new_subline_trimmed = QgsLineString([qgs_pnt_i_trimmed, qgs_pnt_j_trimmed])
            else:
                qgs_ls_new_subline_trimmed = qgs_ls_new_line.clone()
            self.qgs_geom_new_subline_trimmed = QgsGeometry(qgs_ls_new_subline_trimmed)


        return self.qgs_geom_new_subline_trimmed


class RbResults:

    __slots__ = ('in_nbr_features', 'out_nbr_features', 'nbr_bend_reduced', 'nbr_bend_detected', \
                 'qgs_features_out', 'nbr_hole_del', 'nbr_pol_del', 'nbr_pass', 'is_structure_valid', \
                 'lines_log_info')

    def __init__(self):
        self.in_nbr_features = None
        self.out_nbr_features = None
        self.nbr_bend_reduced = 0
        self.nbr_bend_detected = 0
        self.qgs_features_out = None
        self.nbr_hole_del = 0
        self.nbr_pol_del = 0
        self.nbr_pass = 0
        self.is_structure_valid = None
        self.lines_log_info = []

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

class Epsilon():

    __slots__ = '_zero_relative', '_zero_absolute', '_zero_angle', '_map_range'

    def __init__(self, features):
        b_box = features[0].geometry().boundingBox()
        for feature in features:
            b_box.combineExtentWith(feature.geometry().boundingBox())

        delta_x = abs(b_box.xMinimum()) + abs(b_box.xMaximum())
        delta_y = abs(b_box.yMinimum()) + abs(b_box.yMaximum())
        dynamic_xy = max(delta_x, delta_y)
        log_loss = int(math.log(dynamic_xy, 10)+1)
        max_digit = 15  # Number of significative digits for real number
        security = 2
        abs_digit = max_digit  - security
        rel_digit = max_digit - log_loss - security
        self._zero_relative = (1. / (10**(rel_digit)))
        self._zero_absolute = 1. / (10**(abs_digit))
        self._zero_angle= math.radians(.0001)
        self._map_range = max(b_box.width(), b_box.height())

        return

    def set_class_variables(self):

        ReduceBend.ZERO_RELATIVE =  self._zero_relative
        ReduceBend.ZERO_ABSOLUTE = self._zero_absolute
        ReduceBend.ZERO_ANGLE = self._zero_angle
        ReduceBend.MAP_RANGE = self._map_range

        return




#        bends_small = []
#        bends_out = []
#        bends_in = []
#        bends_other = []
#        bends = [(bend.area, i) for i, bend in enumerate(rb_geom.bends) if bend.area < min_adj_area]
#        bends.sort(key=lambda item: item[0])
#        for bend_area, i in bends:
#            if rb_geom.bends[i].area <= min_adj_area * .15:
#                bends_small.append((bend_area,i))
#            elif rb_geom.bends[i].area <= min_adj_area:
#                if rb_geom.bends[i].orientation == 'IN':
#                    bends_in.append((bend_area,i))
#                else:
#                    bends_out.append((bend_area,i))
#            else:
#                bends_other.append((bend_area,i))
#
#        bends = bends_small + bends_out + bends_in + bends_other

#    else:
#
#        bends = [(bend.adj_area, i) for i, bend in enumerate(rb_geom.bends) if bend.area < min_adj_area]
#        bends.sort(key=lambda item: item[0])


class ReduceBend():

    ZERO_ABSOLUTE = None  # To be set later
    ZERO_RELATIVE = None  # To be set later
    ZERO_ANGLE = None     # To be set later
    MAP_RANGE = None      # To be set later

    @staticmethod
    def reduce(qgs_in_features, diameter_tol, nbr_sub_space=1, feedback=None, flag_del_outer=False,
               flag_del_inner=False, validate_structure=False):

        rb = ReduceBend(qgs_in_features, diameter_tol, nbr_sub_space, feedback, flag_del_outer, flag_del_inner,
                        validate_structure)
        results = rb.reduce_bends()

        return results

    @staticmethod
    def extract_polygon_attributes(qgs_geom):
        qgs_line_string = qgs_geom.constGet()
        qgs_pol = QgsPolygon(qgs_line_string.clone())
        area = qgs_pol.area()
        perimeter = qgs_pol.perimeter()

        return (area, perimeter)

    @staticmethod
    def calculate_adj_area(area, perimeter):
        compactness_index = 4 * area * math.pi / perimeter ** 2
        adj_area = area * (.75 / compactness_index)

        return adj_area

    @staticmethod
    def calculate_min_adj_area(diameter_tol):
        min_adj_area = .75 * math.pi * (diameter_tol / 2.) ** 2

        return min_adj_area

    __slots__ = ('qgs_in_features', 'diameter_tol', 'nbr_sub_space', 'feedback', 'flag_del_outer', 'flag_del_inner',
                 'validate_structure', 'rb_collection', 'eps', 'rb_results', 'rb_features', 'rb_geoms')

    def __init__(self, qgs_in_features, diameter_tol, nbr_sub_space, feedback, flag_del_outer, flag_del_inner,
                 validate_structure):

        self.qgs_in_features = qgs_in_features
        self.diameter_tol = diameter_tol
        self.nbr_sub_space = nbr_sub_space
        self.feedback = feedback
        self.flag_del_outer = flag_del_outer
        self.flag_del_inner = flag_del_inner
        self.validate_structure = validate_structure
        self.rb_collection = None

    def reduce_bends(self):
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

#        import cProfile, pstats, io
#        from pstats import SortKey
#        pr = cProfile.Profile()
#        pr.enable()

        self.eps = Epsilon(self.qgs_in_features)
        self.eps.set_class_variables()

        self.rb_results = RbResults()

        # Create the list of RbLineString and RbPoint to process
        self.rb_features = self.create_rb_feature()
        self.rb_results.in_nbr_features = len(self.qgs_in_features)

        # Pre process the LineString: remove to close point and co-linear points
        self.rb_geoms = self.pre_reduction_process()

        # Create the RbCollection a spatial index to accelerate search
        self.rb_collection = RbCollection(self.rb_results, self.nbr_sub_space)
        self.rb_collection.add_features(self.rb_geoms)

        # Execute the bend reduction
        self._manage_reduce_bend()

        # Recreate the QgsFeature
        qgs_features_out = [rb_feature.get_qgs_feature() for rb_feature in self.rb_features]

        # Calculate return values
        self.rb_results.out_nbr_features = len(qgs_features_out)
        self.rb_results.qgs_features_out = qgs_features_out

        # Validate inner spatial structure. For debug purpose only
        if self.validate_structure:
            self.rb_results.is_structure_valid = self.validate_integrity()

#        pr.disable()
#        s = io.StringIO()
#        sortby = SortKey.CUMULATIVE
#        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#        ps.print_stats()
#        print(s.getvalue())

        return self.rb_results

    def create_rb_feature(self):
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

        for qgs_feature in self.qgs_in_features:
            qgs_geom = qgs_feature.geometry()  # extract the Geometry
            feature_type = qgs_geom.wkbType()

            if RbFeature.is_polygon(feature_type):
                rb_features.append(RbPolygon(qgs_feature))
            elif RbFeature.is_line_string(feature_type):
                rb_features.append(RbLineString(qgs_feature))
            elif RbFeature.is_point(feature_type):
                rb_features.append(RbPoint(qgs_feature))
            else:
                raise Exception("Unsupported GeometryType: {}".format(qgs_geom.wkbType()))

        return rb_features

    def pre_reduction_process(self):

        # Delete the outer or inner ring below the diameter tolerance
        if self.flag_del_outer or self.flag_del_inner:
            self.del_outer_inner_ring()

        # Create the list of RbGeom ==> List of geometry to reduce the bend
        rb_geoms = []
        for rb_feature in self.rb_features:
            rb_geoms += rb_feature.get_rb_geom()


        # Remove co-linear and duplicate nodes
        for rb_geom in rb_geoms:
            if rb_geom.qgs_geom.wkbType() == QgsWkbTypes.LineString:
                if not rb_geom.is_simplest:
                    rb_geom.qgs_geom.removeDuplicateNodes(epsilon=ReduceBend.ZERO_RELATIVE)
                    self.delete_co_linear(rb_geom)

        return rb_geoms


    def del_outer_inner_ring(self):
        # Loop over each rb_features
        for i in reversed(range(len(self.rb_features))):
            if isinstance(self.rb_features[i], RbPolygon):  # Only process Polygon
                min_adj_area = ReduceBend.calculate_min_adj_area(self.diameter_tol)
                for j in reversed(range(len(self.rb_features[i].rb_geom))):
                    area, perimeter = ReduceBend.extract_polygon_attributes(self.rb_features[i].rb_geom[j].qgs_geom)
                    adj_area = self.calculate_adj_area(area, perimeter)
                    if j == 0:
                        # Process the exterior ring
                        if self.flag_del_outer:
                            if adj_area < min_adj_area:
                                del self.rb_features[i]  # Delete the rb_feature
                                self.rb_results.nbr_pol_del += 1
                                break
                    else:
                        # Process an interior ring
                        if self.flag_del_inner:
                            if adj_area < min_adj_area:
                                del self.rb_features[i].rb_geom[j]  # Delete the ring
                                self.rb_results.nbr_hole_del += 1

        return

    def set_bend_direction(self, rb_geom):

        if rb_geom.is_closed and len(rb_geom.bends) >= 1:
            first_bend = rb_geom.bends[0]
            qgs_geom_centroid = first_bend.qgs_geom_bend.centroid()
            qgs_geom_pol = QgsGeometry(QgsPolygon(rb_geom.qgs_geom.constGet().clone(), []))
            if qgs_geom_centroid.within(qgs_geom_pol):
                rb_geom.bends[0].direction = 'IN'
            else:
                rb_geom.bends[0].direction = 'OUT'

            for i in range(1, len(rb_geom.bends)):
                if rb_geom.bends[i - 1].direction == 'IN':
                    rb_geom.bends[i].direction = 'OUT'
                else:
                    rb_geom.bends[i].direction = 'IN'


    def pivot_closed_line(self, rb_geom, diameter_tol):

        if rb_geom.is_closed and rb_geom.need_pivot:
            bend_location = None
            bend_area = 0.0
            for bend in rb_geom.bends:
                if bend.area > bend_area:
                    bend_location = bend
                    bend_area = bend.area
                if bend.j - bend.i >= 4:
                    if bend.area >= ReduceBend.calculate_min_adj_area(diameter_tol):
                        bend_location = bend
                        break
            if bend_location is not None:
                # There is bend candidate for a rotation
                qgs_points = rb_geom.qgs_geom.constGet().points()
                new_start_end = (bend_location.j + bend_location.i) // 2
                new_qgs_points = qgs_points[new_start_end:] + qgs_points[1:new_start_end + 1]
                rb_geom.qgs_geom = QgsGeometry(QgsLineString(new_qgs_points))
                if bend.direction == 'OUT':
                    rb_geom.need_pivot = False  # Priority for bend that are going outside

        return


    def _manage_reduce_bend(self):

        rb_geoms_done = []
        nbr_pass = 0
        min_pass = 1
        nbr_geoms = 100.0 / len(self.rb_geoms) if len(self.rb_geoms) >= 1 else 0
        while True:
            current_diameter_tol = self.diameter_tol * min(nbr_pass+1, min_pass)/float(min_pass)
            self.remove_rb_geoms_done(rb_geoms_done)  # Remove feature done to accelerate process
            if nbr_pass == 0:
                self.feedback.setProgress(1) # Set the progress bar
            else:
                self.feedback.setProgress(max(1, int(len(rb_geoms_done)) * nbr_geoms))
            nbr_bend_reduced = 0
            nbr_bend_detected = 0
            for rb_geom in self.rb_geoms:
                if self.feedback.isCanceled(): break
                nbr_bend_detected += self.manage_bend_creation(nbr_pass, rb_geom, self.diameter_tol)
                self.flag_bend_to_reduce(rb_geom, current_diameter_tol)
                nbr_bend_reduced += self.process_bends(rb_geom)

            str = "Iteration: {}; Bends detected: {}; Bend reduced: {}; Tolerance used: {}"\
                  .format(nbr_pass, nbr_bend_detected, nbr_bend_reduced, current_diameter_tol)
            self.rb_results.lines_log_info.append(str)
            if nbr_pass == 0:
                self.rb_results.nbr_bend_detected = nbr_bend_detected

            # While loop breaking condition
            if nbr_pass > min_pass and nbr_bend_reduced == 0:
                break
            nbr_pass += 1

        # Reset the rb_geoms list
        self.rb_geoms += rb_geoms_done
        self.rb_results.nbr_pass = nbr_pass

        return

    def remove_rb_geoms_done(self, rb_geoms_done):

        for i in reversed(range(len(self.rb_geoms))):
            if self.rb_geoms[i].is_simplest:
                rb_geoms_done.append(self.rb_geoms[i])
                del self.rb_geoms[i]

        return


    def delete_co_linear(self, rb_geom):
        # Delete co-linear vertice with angle of near 0 or 180 degrees
#        qgs_line_string = rb_geom.qgs_geom.constGet()
#        qgs_points = qgs_line_string.points()
#        num_points = len(qgs_points)
#        num_points_remaining = num_points
#        if rb_geom.original_geom_type == QgsWkbTypes.LineString:
#            min_remaining_points = 2  # Minimum number of vertice for a LineString
#        else:
#            min_remaining_points = 4  # Minimum number of vertice for a Polygon
        vertex_ids_to_del = []
        angles = rb_geom.get_angles()
        if rb_geom.is_closed and len(angles) >= 1:
            del angles[0]  # Do not process the start/end points (even if co-linear)
        for i, angle in enumerate(angles):
            if abs(angle - math.pi) <= ReduceBend.ZERO_ANGLE or abs(angle) <= ReduceBend.ZERO_ANGLE:
                # Co-linear point or flat angle delete the current point
                vertex_ids_to_del.append(i+1)
#                num_points_remaining -= 1
#        i = num_points - 2
#        p1_x, p1_y = qgs_points[i].x(), qgs_points[i].y()
#        p2_x, p2_y = qgs_points[i + 1].x(), qgs_points[i + 1].y()
#        while i >= 1 and num_points_remaining >= min_remaining_points:
#            p0_x, p0_y = qgs_points[i - 1].x(), qgs_points[i - 1].y()
#            angle = QgsGeometryUtils.angleBetweenThreePoints(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y)
#            if abs(angle - math.pi) <= ReduceBend.ZERO_ANGLE or abs(angle) <= ReduceBend.ZERO_ANGLE:
#                # Co-linear point or flat angle delete the current point
#                vertex_ids_to_del.append(i)
#                num_points_remaining -= 1
#
#            i -= 1
#            p2_x, p2_y = p1_x, p1_y
#            p1_x, p1_y = p0_x, p0_y

        # Delete co-linerar vertex
        for vertex_id_to_del in reversed(vertex_ids_to_del):
            if self.rb_collection is None:
                rb_geom.qgs_geom.deleteVertex(vertex_id_to_del)
            else:
                self.rb_collection.delete_vertex(rb_geom, [vertex_id_to_del])

        if rb_geom.qgs_geom.length() <= ReduceBend.ZERO_RELATIVE:
            # Something wrong.  do not try to simplify the LineString
            rb_geom.is_simplest = True

        return

    def flag_bend_to_reduce(self, rb_geom, diameter_tol):
        # Minimum adjusted area used to find bend to reduce
        min_adj_area = ReduceBend.calculate_min_adj_area(diameter_tol)

        if rb_geom.is_closed and len(rb_geom.bends) >= 3:
            # The closed line start/end point lie on a bend that do not need to be reduced
            del rb_geom.bends[0]  # Remove the first bend
            del rb_geom.bends[-1]  # Remove the last bend

        if rb_geom.is_closed:
            for bend in rb_geom.bends:
                if bend.direction == "OUT": bend.area *= .66
        lst_bends = [(bend.adj_area, i) for i, bend in enumerate(rb_geom.bends) if bend.area < min_adj_area]
        lst_bends.sort(key=lambda item: item[0])

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

        if len(rb_geom.bends) == 0:
            # No more bends to reduce
            rb_geom.is_simplest = True


    def validate_spatial_constraints(self, ind, rb_geom):

        check_constraints = True
        bend = rb_geom.bends[ind]
        qgs_geom_new_subline = bend.get_new_subline(rb_geom)

        # First: check if the bend reduce line string is an OGC simple line
        # We test with a tiny smaller line to ease the testing and false positive error
        if check_constraints:
            if qgs_geom_new_subline.length() >= ReduceBend.ZERO_RELATIVE:
                qgs_geom_new_subline_trimmed = bend.get_new_subline_trimmed(rb_geom)
                qgs_geom_potentials = self.rb_collection.get_line_segments(rb_geom.id,
                                                                      qgs_geom_new_subline_trimmed.boundingBox())
                for qgs_geom_potential in qgs_geom_potentials:
                    if qgs_geom_new_subline_trimmed.disjoint(qgs_geom_potential):
                        #            if not qgs_geom_new_sub_trim_engine.disjoint(qgs_geom_potential.constGet()):
                        # Everything is OK
                        pass
                    else:
                        # The new sub line intersect the line itself. The result would create a non OGC simple line
                        check_constraints = False
                        break
            else:
                qgs_line_string = qgs_geom_new_subline.constGet()
                x = qgs_line_string.startPoint().x()
                y = qgs_line_string.startPoint().y()
                text = "Possibly non OGC simple feature at {},{} use Fix geometries".format(x,y)
                self.feedback.pushInfo(text)


        # Second: check that the new line does not intersect any other line or points
        if check_constraints:
            qgs_rectangle = bend.qgs_geom_bend.boundingBox()
            #        qgs_geom_engine_new_subline = QgsGeometry.createGeometryEngine(qgs_geom_new_subline.constGet())
            qgs_geom_potentials = self.rb_collection.get_features(qgs_rectangle, [rb_geom.id])
            for qgs_geom_potential in qgs_geom_potentials:
                if not qgs_geom_potential.disjoint(qgs_geom_new_subline):
                    #            if not qgs_geom_new_subline.disjoint(qgs_geom_potential):
                    #            if not qgs_geom_engine_new_subline.disjoint(qgs_geom_potential.constGet()):
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

    def process_bends(self, rb_geom):

        nbr_bend_reduced = 0
        for ind in reversed(range(len(rb_geom.bends))):
            bend = rb_geom.bends[ind]
            if bend.to_reduce:
                # Check spatial constraints
                spatial_constraints = self.validate_spatial_constraints(ind, rb_geom)
                if spatial_constraints:
                    v_ids_to_del = list(range(bend.i + 1, bend.j))  # List of the vertex id to delete
                    self.rb_collection.delete_vertex(rb_geom, v_ids_to_del)
                    #                bend.reduce(rb_geom)
                    self.rb_results.nbr_bend_reduced += 1  # Global counter of bend reduced
                    nbr_bend_reduced += 1  # Local counter of bend reduced

        return nbr_bend_reduced

    def create_polygon (self, i, j, qgs_points):

        # Create the list of point to create
        if i < j:
            index = list(range(i, j+1)) + [i]
        else:
            index = list(range(i, len(qgs_points))) + list(range(0,j+1)) + [i]  # Manage circular array

        qgs_sub_points = [qgs_points[k] for k in index]
        qgs_polygon = QgsPolygon(QgsLineString(qgs_sub_points))

        return qgs_polygon


    def manage_bend_creation(self, nbr_pass, rb_geom, diameter_tol):

        self.delete_co_linear(rb_geom)
        nbr_bends_detected = self.detect_bends(rb_geom)
        if rb_geom.is_closed:
            if nbr_pass >= 2 or len(rb_geom.bends) < 10:
                self.pivot_closed_line(rb_geom, diameter_tol)
                self.delete_co_linear(rb_geom)
                nbr_bends_detected = self.detect_bends(rb_geom)

        return nbr_bends_detected

    def detect_bends(self, rb_geom):

        rb_geom.bends = []  # Reset the list of bends
        angles = rb_geom.get_angles()
        # Modify the angle to binary orientation: clockwise or anti clockwise
        orientation = [CLOCK_WISE if angle >= math.pi else ANTI_CLOCK_WISE for angle in angles]
        if rb_geom.is_closed:
            if len(set(orientation)) == 1:
                orientation = []  # All the angles have the same orientation.  No bend to reduce
            else:
                del orientation[0]  # Do not process the fist angle as it is the angle of the start/end

        if len(orientation) >= 1:
            orientation.insert(0, ANTI_CLOCK_WISE) if orientation[0] == CLOCK_WISE else orientation.insert(0, CLOCK_WISE)
            orientation.append(ANTI_CLOCK_WISE) if orientation[-1] == CLOCK_WISE else orientation.append(CLOCK_WISE)

        #            if rb_geom.is_closed:
#                # Copy the first orientation to the last position for the circularity
#                orientation.append(orientation[0])
#            else:
#                orientation.insert(0, ANTI_CLOCK_WISE) if orientation[0] == CLOCK_WISE else orientation.insert(0, CLOCK_WISE)
#                orientation.append(ANTI_CLOCK_WISE) if orientation[-1] == CLOCK_WISE else orientation.append(CLOCK_WISE)

        # Find the inflexion points in the line.
        inflexion = [i for i in range(0, len(orientation)-1) if orientation[i] != orientation[(i + 1)]]
        if len(inflexion) != 0:
#            if rb_geom.is_closed:
#                inflexion.append(inflexion[0])

            qgs_points = rb_geom.qgs_geom.constGet().points()
            for k in range(len(inflexion)-1):
                i = inflexion[k]
                j = (inflexion[(k+1)]+1)
                qgs_polygon = self.create_polygon(i, j, qgs_points)
                rb_geom.bends.append(Bend(i, j, qgs_polygon))

        else:
            # If there is no inflexion the line cannot be simplified
            rb_geom.is_simplest = True

        # Set the direction of the bend inside or outside the polygon
        if rb_geom.is_closed:
            self.set_bend_direction(rb_geom)
        #        if len(rb_geom.bends) == 0:
#            # A line with no inflexion cannot be simplified more
#            rb_geom.is_simplest = True

#        if rb_geom.is_closed:
#            # Process the line string as a closed line
#            angle = QgsGeometryUtils.angleBetweenThreePoints(xy[-1][0],xy[-1][1],xy[0][0],xy[0][1],xy[1][0],xy[1][1])
#            # Add a first orientation linking the start and end of the line (circular array)
#            orientation.insert(0, CLOCK_WISE) if angle >= math.pi else orientation.insert(0, ANTI_CLOCK_WISE)
#            # Find the inflexion points in the line.  Managing the circular array
#            inflexion = [i for i in range(0, len(orientation)) if orientation[i] != orientation[(i+1)%len(orientation)]]
#            for k in range(len(inflexion)):
#                i = inflexion[k]
#                j = inflexion[(k+1)%len(inflexion)]+1
#                if i != j:
#                    qgs_polygon = self.create_polygon(i, j, xy)
#                    rb_geom.bends.append(Bend(i, j, qgs_polygon))
#                else:
#                    # Extemely rare case. Happen with special case polygon
#                    pass

#            if 2 <= len(rb_geom.bends) <= 4:
#                del rb_geom.bends[-1]
#                del rb_geom.bends[-1]
#        else:
#            # For open lines add a point at the start and one at the end (to facilitate )
#            if len(orientation) >= 1:
#                orientation.insert(0, ANTI_CLOCK_WISE) if orientation[0] == CLOCK_WISE else orientation.insert(0, CLOCK_WISE)
#                orientation.append(ANTI_CLOCK_WISE) if orientation[-1] == CLOCK_WISE else orientation.append(CLOCK_WISE)

#            # Find the inflexion points in the line
#            inflexions = [i for i in range(0,len(orientation)-1) if orientation[i] != orientation[i+1]]

#            for k in range(len(inflexions)-1):
#                i = inflexions[k]
#                j = inflexions[k+1]+1
#                qgs_polygon = self.create_polygon(i, j, xy)
#                rb_geom.bends.append(Bend(i, j, qgs_polygon))

        return len(rb_geom.bends)


    def validate_integrity(self):

        is_structure_valid = True
        for rb_geom in self.rb_geoms:
            qgs_line_string = rb_geom.qgs_geom.constGet()
            if qgs_line_string.wkbType() == QgsWkbTypes.LineString:
                qgs_points = qgs_line_string.points()
                for i in range(len(qgs_points)-1):
#                    qgs_rectangle = QgsRectangle(-sys.float_info.max, -sys.float_info.max,
#                                                sys.float_info.max, sys.float_info.max)
#                    spatial_index = self.rb_collection._geom_line_segment._dict_index[5]
#                    b_boxes = spatial_index.intersects(qgs_rectangle)
                    self.rb_collection._geom_line_segment.delete_line_segment(rb_geom.id, qgs_points[i], qgs_points[i+1])

        if is_structure_valid:
            for spatial_index in self.rb_collection._geom_line_segment._dict_index.values():
                qgs_rectangle = QgsRectangle(-sys.float_info.max, -sys.float_info.max,
                                             sys.float_info.max,sys.float_info.max)
                b_boxes = spatial_index.intersects(qgs_rectangle)
#                qgs_pol = spatial_index.geometry(id=b_boxes[0])
                if len(b_boxes) != 0:
                    is_structure_valid = False
                    print ("error: 1")
                    break

        if is_structure_valid:
            for rb_geom in self.rb_geoms:
                qgs_line_string = rb_geom.qgs_geom.constGet()
                qgs_geom = QgsGeometry(qgs_line_string.clone())
                id_b_boxes = self.rb_collection._spatial_index.intersects(qgs_line_string.boundingBox())
                for id_b_box in id_b_boxes:
                    if self.rb_collection._dict_geom_b_boxes [id_b_box] == rb_geom.id:
                        qgs_pol = self.rb_collection._spatial_index.geometry(id_b_box)
                        qgs_geom = qgs_geom.difference(qgs_pol)
                if not qgs_geom.isEmpty():
                    is_structure_valid = False
                    print("error: 2")

        return is_structure_valid
