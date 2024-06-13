import arcpy
import os
import glob

workspace = r'C:\Computation\CUNI\02_LS1\MATKARTO\testing_py_globes\testing\\'
os.chdir(workspace)


# Check if project exists, remove it if it does
if os.path.exists("project.aprx"):
    os.remove("project.aprx")

aprx = arcpy.mp.ArcGISProject(workspace)
# Create project and layout
aprx = arcpy.mp.ArcGISProject("project.aprx")
layout = aprx.createLayout(300, 430, 'MILLIMETER')
layout.name = "Layout"

# Define face positions and rotations
face_positions = [(50, 350), (150, 350), (250, 350), (50, 250), (150, 250), (250, 250),
                  (50, 150), (150, 150), (250, 150), (50, 50), (150, 50), (250, 50)]
face_rotations = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# Define thematic layers
thematic_layers = ['layer1.shp', 'layer2.shp',...]  # list of shapefile paths

# Define scale and projection
scale = 1000000  # in meters
projections = ['+proj=longlat +datum=WGS84 +no_defs',  # list of proj4 strings
               '+proj=longlat +datum=WGS84 +no_defs',
              ...]

# Create map frames and add thematic layers
for i, (position, rotation) in enumerate(zip(face_positions, face_rotations)):
    # Create map frame
    map_frame = layout.createMapFrame(position[0], position[1], 50, 50)
    map_frame.name = f"Face {i+1}"

    # Create map and add thematic layer
    map = aprx.createMap(f"Face {i+1}")
    map_frame.map = map
    layer = arcpy.management.MakeFeatureLayer(thematic_layers[i], f"Layer {i+1}")[0]
    map.addLayer(layer)

    # Set projection and scale
    spatial_ref = arcpy.SpatialReference(text=projections[i])
    map.spatialReference = spatial_ref
    map.scale = scale

    # Set extent to fit thematic layer
    desc = arcpy.Describe(thematic_layers[i])
    extent = desc.extent
    map_frame.camera.setExtent(extent)

# Export layout to PDF
layout.exportToPDF("globe_faces.pdf")