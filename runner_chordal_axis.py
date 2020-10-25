# Code example available here: https://anitagraser.com/2019/03/03/stand-alone-pyqgis-scripts-with-osgeo4w/

import processing

#import lib_geosim

from qgis.analysis import QgsNativeAlgorithms
from qgis.core import QgsApplication, QgsFeatureRequest
from processing.core.Processing import Processing
from processing.tools import dataobjects


# Supply path to qgis install location
QgsApplication.setPrefixPath("/usr/bin/qgis", False)

#profile_folder = 'C:\\Users\\berge\\AppData\\Roaming\\QGIS\\QGIS3\\profiles\\test12'
profile_folder = '.'
# Create a reference to the QgsApplication.  Setting the second argument to False disables the GUI.
app = QgsApplication([], False, profile_folder)

# Load providers
app.initQgis()

Processing.initialize()
#context = dataobjects.createContext()
#context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)



# Call the processing script
try:

    msg = processing.run("script:chordal_axis", { 'INPUT':'/home/daneli/personnel/qgis/geo_file/chordal.gpkg|layername=chord1',
                                                  'CORRECTION':True,
                                                  'OUTPUT':'/home/daneli/personnel/qgis/geo_file/chordal_out1.gpkg'} )
except Exception:
    import traceback
    traceback.print_exc()

#msg = processing.run("script:chordal_axis", {  'INPUT':'C:/Users/berge/PycharmProjects/qgis_geo_sim/road_pol1.gpkg|layername=road_pol1'}, context=context)
#msg = processing.run("script:chordal_axis", { 'INPUT':'C:/Users/berge/PycharmProjects/qgis_geo_sim/road_pol_sub1.gpkg|layername=road_pol_sub'})
#msg = processing.run("script:chordal_axis", { 'INPUT':'C:/Users/berge/PycharmProjects/qgis_geo_sim/test_pol.shp'})
#msg = processing.run("script:chordal_axis", {  'INPUT':'C:/Users/berge/PycharmProjects/qgis_geo_sim/road_pol3.geojson'}, context=context)



print ("End of processing script")

# Stop QGIS appllication
app.exitQgis()
app.exit()