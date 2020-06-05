# Code example available here: https://anitagraser.com/2019/03/03/stand-alone-pyqgis-scripts-with-osgeo4w/

import processing

from qgis.analysis import QgsNativeAlgorithms
from qgis.core import QgsApplication, QgsFeatureRequest
from processing.core.Processing import Processing
from processing.tools import dataobjects

# Supply path to qgis install location
QgsApplication.setPrefixPath("C:\\OSGEO4~1\\apps\\qgis", False)

profile_folder = '.'
# Create a reference to the QgsApplication.  Setting the second argument to False disables the GUI.
app = QgsApplication([], False, profile_folder)

# Load providers
app.initQgis()

Processing.initialize()

# Call the processing script
msg = processing.run("script:sherbend", {'INPUT':r'C:\Users\berge\PycharmProjects\qgis_geo_sim\ttt.gpkg|layername=ttt',
                                         'DIAMETER':5, 'EXCLUDE_HOLE':True, 'EXCLUDE_POLYGON':True,
                                         'OUTPUT':r'C:\Users\berge\PycharmProjects\qgis_geo_sim\road_skeleton_out99.gpkg'} )

print ("End of processing script")

# Stop QGIS appllication
app.exitQgis()
app.exit()