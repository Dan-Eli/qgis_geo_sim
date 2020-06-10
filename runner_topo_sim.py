# Code example available here: https://anitagraser.com/2019/03/03/stand-alone-pyqgis-scripts-with-osgeo4w/

import processing

from qgis.analysis import QgsNativeAlgorithms
from qgis.core import QgsApplication
from processing.core.Processing import Processing

# Supply path to qgis install location
QgsApplication.setPrefixPath("C:\\OSGEO4~1\\apps\\qgis", True)

#profile_folder = 'C:\\Users\\berge\\AppData\\Roaming\\QGIS\\QGIS3\\profiles\\test12'
profile_folder = '.'
# Create a reference to the QgsApplication.  Setting the second argument to False disables the GUI.
app = QgsApplication([], False, profile_folder)

# Load providers
app.initQgis()

Processing.initialize()

# Call the processing script
#msg = processing.run("script:topological_simplifier", {'INPUT':'D:/OneDrive/Personnel/Daniel/QGIS/Kingston/Kingston_dissolved.gpkg|layername=Hydro',
#                                                       'TOLERANCE':3,'TOPOLOGY':True,'OUTPUT':'D:/OneDrive/Personnel/Daniel/QGIS/Kingston/Kingston_dissolved1103.gpkg'})
# Print the name of the algorithm
#QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
#for algo in QgsApplication.processingRegistry().algorithms():
#    print(algo.id(), "------", algo.displayName())

print ("End of processing script")

# Stop QGIS appllication
app.exitQgis()
app.exit()
