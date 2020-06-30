import logging
import sys
from qgis.core import QgsProcessingFeedback

class QgsLoggingFeedback(QgsProcessingFeedback):

    def __init__(self):
        super().__init__()
        self.handler = logging.StreamHandler(sys.stdout)
        # Also show INFO and DEBUG messages
        # self.handler.setLevel(logging.DEBUG)
        logging.getLogger().addHandler(self.handler)

    def reportError(self, msg, fatalError=False):
        logging.log(logging.ERROR, msg)

    def setProgressText(self, text):
        logging.log(logging.INFO, msg)

    def pushInfo(self, info):
        super().pushInfo(info)

    def pushCommandInfo(self, info):
        super().pushCommandInfo(info)

    def pushDebugInfo(self, info):
        super().pushDebugInfo(info)

    def pushConsoleInfo(self, info):
        super().pushConsoleInfo(info)


# Code example available here: https://anitagraser.com/2019/03/03/stand-alone-pyqgis-scripts-with-osgeo4w/

# The import here is not essential but it allows to detect syntax errors
# When not importing these libraries it can be harder to detect error when calling the processing script
import lib_geosim, chordal_axis

import processing
from qgis.analysis import QgsNativeAlgorithms
from qgis.core import QgsApplication, QgsFeatureRequest, QgsLogger
from processing.core.Processing import Processing
from processing.tools import dataobjects






#profile_folder = 'C:\\Users\\berge\\AppData\\Roaming\\QGIS\\QGIS3\\profiles\\test12'
profile_folder = '.'
# Create a reference to the QgsApplication.  Setting the second argument to False disables the GUI.
app = QgsApplication([], False, profile_folder)
# Supply path to qgis install location
QgsApplication.setPrefixPath("C:\\OSGEO4~1\\apps\\qgis", True)
# Load providers
app.initQgis()

Processing.initialize()
#context = dataobjects.createContext()
#context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)



#QgsLogger.debug('This is a debug')
#QgsLogger.warning('This is a warning')
#QgsLogger.critical('This is critical')
#myLogger = QgsLogger("coco.log")
#myLogger.debug('This is a debug')
#myLogger.warning('This is a warning')
#myLogger.critical('This is critical')

from qgis.core import QgsMessageLog
QgsMessageLog.logMessage("File all pixels set to no data" )


# Call the processing script
try:
    feedback = QgsLoggingFeedback()
#    msg = processing.run("script:chordal_axis", { 'INPUT':r'C:\Users\berge\PycharmProjects\qgis_geo_sim\road_pol.gpkg|layername=road_pol',
#                                                  'CORRECTION':False,
#                                                  'OUTPUT':r'C:\Users\berge\PycharmProjects\qgis_geo_sim\road_skeleton_out12326.gpkg'},feedback=feedback )
    msg = processing.run("script:chordal_axis", { 'INPUT':r'C:\Users\berge\PycharmProjects\qgis_geo_sim\small_test.gpkg|layername=small_test',
                                                  'CORRECTION':False,
                                                  'OUTPUT':r'C:\Users\berge\PycharmProjects\qgis_geo_sim\road_skeleton_out12329.gpkg'},feedback=feedback )
except Exception:
    import traceback
    traceback.print_exc()

#msg = processing.run("script:chordal_axis", {  'INPUT':'C:/Users/berge/PycharmProjects/qgis_geo_sim/road_pol1.gpkg|layername=road_pol1'}, context=context)
#msg = processing.run("script:chordal_axis", { 'INPUT':'C:/Users/berge/PycharmProjects/qgis_geo_sim/road_pol_sub1.gpkg|layername=road_pol_sub'})
#msg = processing.run("script:chordal_axis", { 'INPUT':'C:/Users/berge/PycharmProjects/qgis_geo_sim/test_pol.shp'})
#msg = processing.run("script:chordal_axis", {  'INPUT':'C:/Users/berge/PycharmProjects/qgis_geo_sim/road_pol3.geojson'}, context=context)

# Stop QGIS appllication
app.exitQgis()
app.exit()

print ("End of processing script")