#Remaining modifications:
#  - put comment in the code
#  - Manage case where a closed line string is of the same orientation except the start/end which is different
#  - Add some comprehension list instead of for loop
#  - edit line with a smooth line instead of a straight line
#  - reorient the bend solely on the second pass as it will mobe start/end probably on a greater bend
#  - in the bend flagging process prioritize the bend that goes outside the polygon first (for polygon)
#  - edit closed line in comments
#  - put some code to correct J-Bend
#  - ajouter une class Reduce et placer la majorité des routines orphelines
#  - Test reduce_bend performance with the profiler

from abc import ABC, abstractmethod
import math
from qgis.core import QgsFeature, QgsPoint, QgsPointXY, QgsLineString, QgsPolygon, QgsWkbTypes, QgsSpatialIndex, \
    QgsGeometry, QgsGeometryUtils, QgsRectangle

# Define global constant
GeoSim_NBR_PARTITION = 100
GeoSim_EPSILON_REL = None
GeoSim_EPSILON_ABS = None
GeoSim_CW = 0
GeoSim_ACW = -1


def pre_reduction_process(rb_features, rb_results, diameter_tol, flag_del_outer, flag_del_inner):

    # Delete the outer or inner ring below the diameter tolerance
    if flag_del_outer or flag_del_inner:
        del_outer_inner_ring(rb_features, rb_results, diameter_tol, flag_del_outer, flag_del_inner)

    # Create the list of RbGeom ==> List of geometry to reduce the bend
    rb_geoms = []
    for rb_feature in rb_features:
        rb_geoms += rb_feature.get_rb_geom()

    # Remove near duplicate points for each geometry to reduce
    for rb_geom in rb_geoms:
        if rb_geom.qgs_geom.wkbType() == QgsWkbTypes.LineString:
            if rb_geom.qgs_geom.length() > GeoSim_EPSILON_REL:
                rb_geom.qgs_geom.removeDuplicateNodes(epsilon=GeoSim_EPSILON_REL)
            else:
                rb_geom.is_simplest = True  # Something wrong here...

    # Remove co-linear and almost co-linear points
    for rb_geom in rb_geoms:
        if rb_geom.qgs_geom.wkbType() == QgsWkbTypes.LineString:
            if not rb_geom.is_simplest:
                delete_co_linear(None, rb_geom)

    # Set is_simplest attribute according to some basic rules
    for rb_geom in rb_geoms:
        if not rb_geom.is_simplest:
            if rb_geom.original_geom_type == QgsWkbTypes.Point:
                rb_geom.is_simplest = True  #  A point is always at its simplest form
            elif rb_geom.original_geom_type == QgsWkbTypes.LineString:
                if rb_geom.qgs_geom.constGet().numPoints() <= 2:  # A LineString must have 3 or more points to be candidate for reduction
                    rb_geom.is_simplest = True
            elif rb_geom.original_geom_type == QgsWkbTypes.Polygon:
                if rb_geom.qgs_geom.constGet().numPoints() <= 4:  # A Polygon must have 3 or more points to be candidate for reduction
                    rb_geom.is_simplest = True
            else:
                raise Exception ("Unhandled geometry type: {}".format(QgsWkbTypes.displayString(rb_geom.original_geom_type)))

    return rb_geoms


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

    eps = Epsilon(qgs_in_features)

    rb_results = RbResults()

    # Create the list of RbLineString and RbPoint to process
    rb_features = create_rb_feature(qgs_in_features)
    rb_results.in_nbr_features = len(qgs_in_features)

    # Pre process the LineString: remove to close point and co-linear points
    rb_geoms = pre_reduction_process(rb_features, rb_results, diameter_tol, flag_del_outer, flag_del_inner)


    # Create the RbCollection a spatial index to accelerate search
    rb_collection = RbCollection()
    rb_collection.add_features(rb_geoms)

    # Execute the bend reduction
    _manage_reduce_bend(rb_geoms, rb_collection, rb_results, diameter_tol, feedback)

    # Recreate the QgsFeature
    qgs_features_out = [rb_feature.get_qgs_feature() for rb_feature in rb_features]

    # Calculate return values
    rb_results.out_nbr_features = len(qgs_features_out)
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
            qgs_features = []
            for i in range(len(qgs_points)-1):
                qgs_feature= QgsFeature(id=self._get_next_id())
                qgs_feature.setGeometry(QgsLineString(qgs_points[i], qgs_points[i+1]))
                qgs_features.append(qgs_feature)

            # Add the line segment in the spatial index
            spatial_index.addFeatures(qgs_features)
            self._dict_index[rb_geom.id] = spatial_index

        return

    def delete_line_segment(self, qgs_geom_id, qgs_pnt0, qgs_pnt1):

        spatial_index = self._dict_index[qgs_geom_id]
#        qgs_rectangle = QgsRectangle(QgsPointXY(qgs_pnt0), QgsPointXY(qgs_pnt1))
        qgs_mid_point = QgsGeometryUtils.midpoint(qgs_pnt0, qgs_pnt1)
        qgs_rectangle = qgs_mid_point.boundingBox()
        qgs_rectangle.grow(GeoSim_EPSILON_REL*100)
#        qgs_rectangle.grow(100*GeoSim_EPSILON_REL)
        ids = spatial_index.intersects(qgs_rectangle)
        for id in ids:
            qgs_geom_line = spatial_index.geometry(id)  # Extract geometry from index
            qgs_pnt_start = qgs_geom_line.vertexAt(0)
            qgs_pnt_end = qgs_geom_line.vertexAt(1)
            if qgs_pnt_start.distance(qgs_pnt0) <= GeoSim_EPSILON_REL and \
               qgs_pnt_end.distance(qgs_pnt1) <= GeoSim_EPSILON_REL:
                feature = QgsFeature(id=id)
                feature.setGeometry(QgsLineString([qgs_pnt_start, qgs_pnt_end]))
                if (spatial_index.deleteFeature(feature)):
                    deleted = True
                    break
                else:
                    0/0
            else:
                deleted = False

        if not deleted:
            0/0

        return

    def add_line_segment(self, rb_geom_id, qgs_pnt0, qgs_pnt1):

#        qgs_line_segment = QgsLineString(qgs_pnt0, qgs_pnt1)
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

    def __init__(self):
        self._spatial_index = QgsSpatialIndex(flags=QgsSpatialIndex.FlagStoreFeatureGeometries)
        self._dict_qgs_geom = {}
        self._dict_geom_b_boxes = {}
        self._id_b_box = 0
        self._geom_line_segment = GeomLineSegment()

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

        if max(qgs_b_box.width(),qgs_b_box.height()) <= self.min_cell_size:
            # No bounding box should be 0 or near 0 in width or in height
            # Also avoid some side effects when bends are parallel to the boundind box
            qgs_b_box.grow(GeoSim_EPSILON_REL*100)
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
            if x_max-x_min < self.min_cell_size:
                # Only create north-south quadrant
                qgs_rect_s = QgsRectangle(x_min, y_min, x_max, y_mid)
                qgs_rect_n = QgsRectangle(x_min, y_mid, x_max, y_max)
                qgs_rectangles = [qgs_rect_s, qgs_rect_n]
            elif y_max-y_min < self.min_cell_size:
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

        self.min_cell_size = GeoSim_MAP_RANGE/GeoSim_NBR_PARTITION
        dummy = 0
        for rb_geom in rb_geoms:
            dummy += self.create_index_b_boxes(rb_geom.id, rb_geom.qgs_geom)
#            qgs_features_b_box = self._subdivide_b_box(rb_geom.qgs_geom, rb_geom.qgs_geom.boundingBox(), [])
#            b_box_ids = [qgs_feature_b_box.id() for qgs_feature_b_box in qgs_features_b_box]
#            for b_box_id in b_box_ids:
#                dummy += len(b_box_ids)
#                self._dict_geom_b_boxes[b_box_id] = rb_geom.id  # Keep link between b_box id and the rb_geom id

#            self._spatial_index.addFeatures(qgs_features_b_box)  # Add the bounding boxes into the index
            self._dict_qgs_geom[rb_geom.id] = rb_geom.qgs_geom  # Keep a reference to the qgs_geom
            self._geom_line_segment.add_geom(rb_geom)  #  Create the line segment index for the geometry

        print ("Nbr bbox: ", dummy)

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

#            print("need to add the line segment bounding box in the spatial index")

        # Now that the line segment spatial structure is updated let's delete the vertex
        for v_id_to_del in reversed(v_ids_to_del):
            qgs_geom.deleteVertex(v_id_to_del)

        return


    def get_line_segments(self, geom_id, qgs_rectangle):

        qgs_geom_line_segments = self._geom_line_segment.get_line_segments(geom_id, qgs_rectangle)

        return qgs_geom_line_segments


class RbGeom:
    _id_counter = 0

    def __init__(self, qgs_abs_geom, original_geom_type):

        self.id = self.next_id()
        self.original_geom_type = original_geom_type
        self.is_simplest = False
        self.qgs_geom = QgsGeometry(qgs_abs_geom)
        self.qgs_rectangle = self.qgs_geom.boundingBox()
        self.bends = []
        self.nbr_bend_reduced = 0


    def edit_closed_line(self, diameter_tol):

        0/0
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

    def next_id(self):

        RbGeom._id_counter += 1

        return (RbGeom._id_counter)


class Bend:

    def __init__(self, i, j, qgs_points):

        self.i = i
        self.j = j
        self.qgs_geom_bend = QgsGeometry(QgsPolygon(QgsLineString(qgs_points[i:j + 1])))
        self.area = self.qgs_geom_bend.area()
        self.perimeter = self.qgs_geom_bend.length()
###        self.qgs_geom_bend = QgsGeometry(qgs_pol_bend.clone())
        self.adj_area = calculate_adj_area(self.area, self.perimeter)
        self.to_reduce = False
        self.qgs_geom_new_subline = None
        self.qgs_geom_new_subline_trimmed = None

#    def reduce(self, rb_geom):
#
#        for ind in range(self.j - 1, self.i, -1):  # Process the vertice in reverse order
#            rb_geom.qgs_geom.deleteVertex(ind)  # Delete vertex to reduce line
#
#        # Update statistics
#        rb_geom.nbr_bend_reduced += 1

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
            qgs_pnt_i_trimmed = qgs_ls_new_line.interpolatePoint(GeoSim_EPSILON_REL)
            qgs_pnt_j_trimmed = qgs_ls_new_line.interpolatePoint(qgs_ls_new_line.length() - GeoSim_EPSILON_REL)
            qgs_ls_new_subline_trimmed = QgsLineString([qgs_pnt_i_trimmed, qgs_pnt_j_trimmed])
            self.qgs_geom_new_subline_trimmed = QgsGeometry(qgs_ls_new_subline_trimmed)

        return self.qgs_geom_new_subline_trimmed


class RbResults:

    def __init__(self):
        self.in_nbr_features = None
        self.out_nbr_features = None
        self.nbr_bend_reduced = 0
        self.nbr_bend_detected = 0
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

class Epsilon():

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
        abs_digit = max_digit  # Number of significative number in real*8
        rel_digit = max_digit - log_loss - security
        global GeoSim_EPSILON_REL, GeoSim_EPSILON_ABS, GeoSim_MAP_RANGE
        GeoSim_EPSILON_REL = (1. / (10**(rel_digit)))
        GeoSim_EPSILON_ABS = 1. / (10**(abs_digit))
        GeoSim_MAP_RANGE = max(b_box.width(), b_box.height())

        return


def create_rb_feature(qgs_features):
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

        if RbFeature.is_polygon(feature_type):
            rb_features.append(RbPolygon(qgs_feature))
        elif RbFeature.is_line_string(feature_type):
            rb_features.append(RbLineString(qgs_feature))
        elif RbFeature.is_point(feature_type):
            rb_features.append(RbPoint(qgs_feature))
        else:
            raise Exception("Unsupported GeometryType: {}".format(qgs_geom.wkbType()))

    return rb_features


#def angle_between_three_points(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y):
#
#    angle1 = math.atan2(p0_y - p1_y, p0_x - p1_x )
#    angle2 = math.atan2(p2_y - p1_y, p2_x - p1_x)
#
#    angle = angle1 - angle2
#
#    # Normalizing angle
#    clippedAngle = angle
#    if clippedAngle >= math.pi * 2 or clippedAngle <= -2 * math.pi:
#        clippedAngle = clippedAngle % 2*math.pi
#    if clippedAngle < 0.0:
#        clippedAngle += 2 * math.pi
#
#    return clippedAngle

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


def delete_co_linear(rb_collection, rb_geom):
    # Delete almost duplicate point
#    qgs_line_string = rb_geom.qgs_geom.constGet()
#    qgs_line_string.removeDuplicateNodes(GeoSim_EPSILON)

    # Delete co-linear points on the line angle of near 0 or 180 degrees
#    qgs_points = qgs_line_string.points()
#    coords_xy = create_coords_xy(rb_geom.qgs_geom)
    qgs_line_string = rb_geom.qgs_geom.constGet()
    qgs_points = qgs_line_string.points()
#    coords_xy = [(qgs_point.x(),qgs_point.y()) for qgs_point in qgs_line_string.points()]
    num_points = len(qgs_points)
    num_points_remaining = num_points
    if rb_geom.original_geom_type == QgsWkbTypes.LineString:
        min_remaining_points = 2  # Minimum number of vertice for a LineString
    else:
        min_remaining_points = 4 # Minimum number of vertice for a Polygon
#    num_points = len(coords_xy)
    vertex_ids_to_del = []
    i = num_points - 2
    p1_x, p1_y = qgs_points[i].x(), qgs_points[i].y()
    p2_x, p2_y = qgs_points[i + 1].x(), qgs_points[i + 1].y()
    while i >= 1 and num_points_remaining >= min_remaining_points:
        p0_x, p0_y = qgs_points[i-1].x(), qgs_points[i-1].y()
#        p0_x, p0_y = coords_xy[i - 1][0], coords_xy[i - 1][1]
#        p1_x, p1_y = coords_xy[i][0], coords_xy[i][1]
#        p2_x, p2_y = coords_xy[i + 1][0], coords_xy[i + 1][1]

        angle = QgsGeometryUtils.angleBetweenThreePoints(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y)
#        angle = angle_between_three_points(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y)
        #        if abs(angle-angle_prime) > GeoSim_EPSILON:
        #            print ("calcul angle pas bon...")
        #            0/0
        if abs(angle - math.pi) <= GeoSim_EPSILON_ABS or abs(angle) <= GeoSim_EPSILON_ABS:
            # Co-linear point or flat angle delete the current point
            #del coords_xy[i]
#            del coords_xy[i]
            vertex_ids_to_del.append(i)
            num_points_remaining -= 1

        i -= 1
        p2_x, p2_y = p1_x, p1_y
        p1_x, p1_y = p0_x, p0_y

    # Delete co-linerar vertex
    for vertex_id_to_del in vertex_ids_to_del:
        if rb_collection is None:
            rb_geom.qgs_geom.deleteVertex(vertex_id_to_del)
        else:
            rb_collection.delete_vertex(rb_geom, [vertex_id_to_del])


    if rb_geom.qgs_geom.length() <= GeoSim_EPSILON_REL:
        # Something wrong.  do not try to simplify the LineString
        rb_geom.is_simplest = True

    return


def detect_bends(rb_geom):

#    # Delete co-linear and almost duplicate point
#    delete_co_linear(rb_collection, rb_geom)

    qgs_line_string = rb_geom.qgs_geom.constGet()
    qgs_points = qgs_line_string.points()
#    coords_xy = [(qgs_point.x(), qgs_point.y()) for qgs_point in qgs_points]
#    if len(qgs_points) != len(coords_xy):
#        raise Exception ("Internal corruption detected in module: detect_bends")

#    coords_xy = [(qgs_point.x(), qgs_point.y()) for qgs_point in qgs_line_string.points()]

    num_qgs_points = len(qgs_points)
    i = 1
    angles = []
    p0_x, p0_y = qgs_points[i-1].x(), qgs_points[i-1].y()
    p1_x, p1_y = qgs_points[i].x(), qgs_points[i].y()
    # Extract the angles at each vertice except start and end vertice
    while i <= num_qgs_points - 2 and num_qgs_points >= 3:
        p2_x, p2_y = qgs_points[i+1].x(), qgs_points[i+1].y()

#        p0_x, p0_y = coords_xy[i - 1][0], coords_xy[i - 1][1]
#        p1_x, p1_y = coords_xy[i][0], coords_xy[i][1]
#        p2_x, p2_y = coords_xy[i + 1][0], coords_xy[i + 1][1]
        angles.append(QgsGeometryUtils.angleBetweenThreePoints(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y))
        #        angle_prime = angle_between_three_points(p0_x, p0_y, p1_x, p1_y, p2_x, p2_y)
        #        if abs(angles[-1] - angle_prime) > GeoSim_EPSILON:
        #            print("calcul angle pas bon...")
        #            0 / 0
        p0_x, p0_y = p1_x, p1_y
        p1_x, p1_y = p2_x, p2_y
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

        if qgs_line_string.isClosed() and i + 2 - start == num_qgs_points:
            # A closed lines with only one bend cannot be reduced
            rb_geom.bends = []  # Reset the bends

    if len(rb_geom.bends) == 0:
        # There is no bend so nothing to simplify
        rb_geom.is_simplest = True

#    if rb_geom.nbr_bend_detected is None:
#        # For the first pass set the number of bend detected (for statistics purpose)
#        rb_geom.nbr_bend_detected = len(rb_geom.bends)

    return len(rb_geom.bends)


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


def priorize_bend_reduction(rb_geom, min_adj_area):

    if rb_geom.qgs_geom.constGet().isClosed() and len(rb_geom.bends) >=1:
        first_bend = rb_geom.bends[0]
        qgs_geom_centroid = first_bend.qgs_geom_bend.centroid()
        qgs_geom_pol = QgsGeometry(QgsPolygon(rb_geom.qgs_geom.constGet().clone(), []))
        if qgs_geom_centroid.within(qgs_geom_pol):
            rb_geom.bends[0].orientation = 'IN'
        else:
            rb_geom.bends[0].orientation= 'OUT'

        for i in range(1,len(rb_geom.bends)):
            if rb_geom.bends[i-1].orientation == 'IN':
                rb_geom.bends[i].orientation = 'OUT'
            else:
                rb_geom.bends[i].orientation = 'IN'

        bends_small = []
        bends_out = []
        bends_in = []
        bends_other = []
        bends = [(bend.area, i) for i, bend in enumerate(rb_geom.bends) if bend.area < min_adj_area]
        bends.sort(key=lambda item: item[0])
        for bend_area, i in bends:
            if rb_geom.bends[i].area <= min_adj_area * .15:
                bends_small.append((bend_area,i))
            elif rb_geom.bends[i].area <= min_adj_area:
                if rb_geom.bends[i].orientation == 'IN':
                    bends_in.append((bend_area,i))
                else:
                    bends_out.append((bend_area,i))
            else:
                bends_other.append((bend_area,i))

        bends = bends_small + bends_out + bends_in + bends_other

    else:

        bends = [(bend.adj_area, i) for i, bend in enumerate(rb_geom.bends) if bend.area < min_adj_area]
        bends.sort(key=lambda item: item[0])

    return bends


    return bends


def flag_bend_to_reduce(rb_geom, diameter_tol):
    # Minimum adjusted area used to find bend to reduce
    min_adj_area = calculate_min_adj_area(diameter_tol)

    if rb_geom.qgs_geom.constGet().isClosed() and len(rb_geom.bends) >= 3:
        # The closed line has been rotated and the start/end point lie on a bend that do not need to be reduced
        del rb_geom.bends[0]  # Remove the first bend
        del rb_geom.bends[-1]  # Remove the last bend

#    lst_bends = priorize_bend_reduction(rb_geom, min_adj_area)
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

#    if len(lst_bends) == 0 or lst_bends[0][0] >= min_adj_area:
#        rb_geom.is_simplest = True


def validate_spatial_constraints(ind, rb_geom, rb_collection):

    check_constraints = True
    bend = rb_geom.bends[ind]
#    qgs_geom_line_string = rb_geom.qgs_geom

    # First: check if the bend reduce line string is an OGC simple line
    # We test with a tiny smaller line to ease the testing and false positive error
    if check_constraints:
        qgs_geom_new_subline_trimmed = bend.get_new_subline_trimmed(rb_geom)
#        qgs_geom_new_sub_trim_engine = QgsGeometry.createGeometryEngine(qgs_geom_new_subline_trimmed.constGet())
        qgs_geom_potentials = rb_collection.get_line_segments(rb_geom.id, qgs_geom_new_subline_trimmed.boundingBox())
        for qgs_geom_potential in qgs_geom_potentials:
            if qgs_geom_new_subline_trimmed.disjoint(qgs_geom_potential):
#            if not qgs_geom_new_sub_trim_engine.disjoint(qgs_geom_potential.constGet()):
                # Everything is OK
                pass
            else:
                # The new sub line intersect the line itself. The result would create a non OGC simple line
                check_constraints = False
                break

    # Second: check that the new line does not intersect any other line or points
    if check_constraints:
        qgs_rectangle = bend.qgs_geom_bend.boundingBox()
        qgs_geom_new_subline = bend.get_new_subline(rb_geom)
#        qgs_geom_engine_new_subline = QgsGeometry.createGeometryEngine(qgs_geom_new_subline.constGet())
        qgs_geom_potentials = rb_collection.get_features(qgs_rectangle, [rb_geom.id])
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


def process_bends(rb_collection, rb_results, rb_geom):

#    for pass_num in (0,1):
#        nbr_bend_reduced = 0
#        if pass_num == 0:
#            if rb_geom.qgs_geom.isSimple():
#                qgs_points_copy = rb_geom.qgs_geom.constGet().points()
#                validate_is_simple = False
#            else:
#                # The line is already not simple.  Should not happen but algothim managed this case
#                validate_is_simple = True
#                print("OGC non valide avant pass 0")
#        else:
#            # A second pass indicates a problem of line self intersecting (manage intersection of each bend reduction)
#            validate_is_simple = True

    nbr_bend_reduced = 0
    for ind in reversed(range(len(rb_geom.bends))):
        bend = rb_geom.bends[ind]
        if bend.to_reduce:
            # Check spatial constraints
            spatial_constraints = validate_spatial_constraints(ind, rb_geom, rb_collection)
            if spatial_constraints:
                v_ids_to_del = list(range(bend.i+1,bend.j))  # List of the vertex id to delete
                rb_collection.delete_vertex(rb_geom, v_ids_to_del)
#                bend.reduce(rb_geom)
                rb_results.nbr_bend_reduced += 1 # Global counter of bend reduced
                nbr_bend_reduced += 1  # Local counter of bend reduced

#        if pass_num == 0:
#            # Validate for line self intersection
#            if rb_geom.qgs_geom.isSimple():
#                # The second pass is not necessary
#                break
#            else:
#                # Undo the edit and recreate the original QgsLineString
#                qgs_line_string = QgsLineString(qgs_points_copy)
#                rb_geom.qgs_geom = QgsGeometry(qgs_line_string.clone())
#                print("OGC non valide apres pass 0")
#        else:
#            # Nothing to validate in a second pass
#            pass

    return nbr_bend_reduced



def _manage_reduce_bend(rb_geoms, rb_collection, rb_results, diameter_tol, feedback):

    rb_geoms_done = []
    nbr_pass = 0
    previous_pass_nbr_bends = -1
    current_pass_nbr_bends = 0
    nbr_geoms = 100.0 / len(rb_geoms) if len(rb_geoms) >= 1 else 0
    while nbr_pass < 6 or previous_pass_nbr_bends != current_pass_nbr_bends:
        remove_rb_geoms_done(rb_geoms, rb_geoms_done)  # Remove feature done to accelerate process
        # set the progress bar
        if feedback is not None:
            if len(rb_geoms_done) == 0:
                feedback.setProgress(1)
            else:
                feedback.setProgress(int(len(rb_geoms_done) * nbr_geoms))
        previous_pass_nbr_bends = current_pass_nbr_bends
        current_pass_nbr_bends = 0
        if nbr_pass <= 1:
            current_diameter_tol = diameter_tol * .15
        elif nbr_pass <= 3:
            current_diameter_tol = diameter_tol * .5
        else:
            current_diameter_tol = diameter_tol
        for rb_geom in rb_geoms:
            delete_co_linear(rb_collection, rb_geom)
            nbr_bend_detected = detect_bends(rb_geom)
            if nbr_pass == 0:
                rb_results.nbr_bend_detected += nbr_bend_detected
#            if nbr_pass == 0:
#                # Edit the start/end point for closed QgsLineString
#                rb_geom.edit_closed_line(diameter_tol)
            flag_bend_to_reduce(rb_geom, current_diameter_tol)
            current_pass_nbr_bends += process_bends(rb_collection, rb_results, rb_geom)
        # Check if all bend are processed
#        is_terminated = is_bend_reduction_terminated(last_nbr_bend_reduced, rb_geoms)
        print ("Passe: ", nbr_pass, "   nbr bends: ", current_pass_nbr_bends)
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