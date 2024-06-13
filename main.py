import arcpy
import numpy as np
import os



# Path to the continents shapefile
continents = r"data\continents.shp"





# Check if the geodatabase exists
if arcpy.Exists(gdb_path):
    # Delete the existing geodatabase
    arcpy.Delete_management(gdb_path)
    print("Existing geodatabase deleted.")
# Create a new file geodatabase
arcpy.CreateFileGDB_management(workspace, gdb_name)
print("New geodatabase created.")


# Create the file geodatabase - Name of the geodatabase
gdb_name = "py_globe_base.gdb"
gdb_path = os.path.join(workspace, gdb_name)


