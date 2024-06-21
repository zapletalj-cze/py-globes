import arcpy
import numpy as np
import os
import glob
import tools
import face_definitions
import sys
from shapely.geometry import Polygon, Point, LineString
from math import pi, sin, cos, tan
from collections import namedtuple


workspace = r"C:\Computation\CUNI\02_LS1\MATKARTO\testing_py_globes"

# DO NOT CHANGE
arcpy.env.workspace = workspace
arcpy.env.overwriteOutput = 1
os.chdir(workspace)
# create processing folder in workspace
if not os.path.exists("processing"):
    os.mkdir("processing")
if not os.path.exists("check"):
    os.mkdir("check")

R = 6378000

# set basemap
# Mapnik_OSM, ESRI_Imagery, Esri_Shaded_Relief, Railway_Map
basemap = 'Mapnik_OSM'

basemaps = {
            'Mapnik_OSM': 'https://tile.openstreetmap.org/{z}/{x}/{y}.png',
            'ESRI_Imagery': 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
            'Esri_Shaded_Relief': 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}',
            'Railway_Map': 'http://sgx.geodatenzentrum.de/wmts_topplus_open/tile/1.0.0/web/default/WEBMERCATOR/{z}/{y}/{x}.png',
            'Google Sattelite': 'https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}',
            'Mapzen Global Terrain': 'https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{z}/{x}/{y}.png'
            }

# load data in shapefile or geojson
data_shp = glob.glob(r"data/*.shp")
data_json = glob.glob(r"data/*.geojson")
data = data_shp + data_json
print(f'Vector datasets found: {data}')
side_len = 50
cir_rad = (side_len / 2) * 1.618

# Open ArcGIS Pro project
project_path = r"project\globes.aprx"
pro_project = arcpy.mp.ArcGISProject(project_path) # Creates a blank project

# Remove all layouts from project
for layout in pro_project.listLayouts():
    layout.removeElement(layout)

# remove all maps from project
maps = pro_project.listMaps()  # Get a list of maps
for map in maps: 
    pro_project.deleteItem(map)  
# create A3 layout
layout = pro_project.createLayout(594, 420, 'MILLIMETER')


# pentagon shifters
def shift_s(d):
    return sin(pi / 5) * d / tan(pi / 5)

def shift_s2(d):
    return sin(2*pi / 5) * d / tan(pi / 5)

def shift_c(d):
    return cos(pi / 5) * d / tan(pi / 5)

def shift_c2(d):
    return cos(2*pi / 5) * d / tan(pi / 5)

def shift_d(d):
    return 2*(tan(0.3*pi) * (d / 2))

def cir_rad(d):
    R = d / (2 * sin(pi / 5))
    return R

# load faces
faces = face_definitions.face_data
# definision for northern hemisphere
shift_x = [0, 
            shift_s2(side_len),
            shift_s(side_len),
            -shift_s(side_len),
            -shift_s2(side_len),
            0, 
            0, 
            0, 
            0, 
            0, 
            0, 
            0]

shift_y = [
            -shift_d(side_len), 
            -shift_c2(side_len), 
            shift_c(side_len), 
            shift_c(side_len), 
            -shift_c2(side_len), 
            0,
            0, 
            0, 
            0, 
            0, 
            0, 
            0]

rotate = [0, 72, 144, -144, -72, 0, -72, -144, 144, 72, 0, 36]

# center for N Africa
val = faces['Face 1']
pentagram = tools.layout.pent_create(a = side_len, angle = val[10])
pentagram_shift = tools.layout.pent_move(pentagram, x_shift=val[8]+shift_x[0], y_shift=val[9]+shift_y[0]) # face 1
ca_c = tools.layout.pent_center(pentagram_shift)
ca_c_x, ca_c_y = ca_c[0], ca_c[1]

print(f'Center_x : {ca_c_x}, center_y: {ca_c_y}')
# update southern hemisphere
Shift = namedtuple("Shift", ["x_value", "x_index", "y_value", "y_index"])
shift_south = [
    Shift(ca_c_x + shift_s(side_len), 5, ca_c_y - shift_c(side_len), 5),
    Shift(ca_c_x + shift_s(side_len), 11, ca_c_y - shift_c(side_len)-shift_d(side_len), 11)
]
for shift in shift_south:
    if 0 <= shift.x_index < len(shift_x):
        shift_x[shift.x_index] = shift.x_value
    if 0 <= shift.y_index < len(shift_y):
        shift_y[shift.y_index] = shift.y_value

# Set center to antlantida
val = faces['Face 12']
pentagram = tools.layout.pent_create(a = side_len, angle = val[10])
pentagram_shift = tools.layout.pent_move(pentagram, x_shift=val[8]+shift_x[11], y_shift=val[9]+shift_y[11]) # face 1
ca_c = tools.layout.pent_center(pentagram_shift)    
ca_c_x, ca_c_y = ca_c[0], ca_c[1]
# shift for southern hemisphere update
Shift = namedtuple("Shift", ["x_value", "x_index", "y_value", "y_index"])
shift_south = [
    Shift(ca_c_x - shift_s(side_len), 8, ca_c_y - shift_c(side_len), 8),
    Shift(ca_c_x + shift_s(side_len), 7, ca_c_y - shift_c(side_len), 7),
    Shift(ca_c_x + shift_s2(side_len), 6, ca_c_y + shift_c2(side_len), 6),
    Shift(ca_c_x - shift_s2(side_len), 9, ca_c_y + shift_c2(side_len), 9),


]

for shift in shift_south:
    if 0 <= shift.x_index < len(shift_x):
        shift_x[shift.x_index] = shift.x_value
    if 0 <= shift.y_index < len(shift_y):
        shift_y[shift.y_index] = shift.y_value


i = 0
# iterrate over keys and values in faces
for k, v in faces.items():
    if i == 2 or i == 10 or i == 0 or i ==4 or i ==1 or i ==3 or i == 5 or i == 11 or i == 9 or i ==8 or i ==7 or i ==6:
        print(f'WORKING ON :{k}')
        # Create new map to project
        map = pro_project.createMap(f"Map_{i+1}")
        for lyr in map.listLayers(): 
            map.removeLayer(lyr)

        ## BOUNDARY
        # create face's boundary
        # convert degrees to radians
        U_B = np.radians(v[0])
        V_B = np.radians(v[1])
        # uk = tools.projections.normalize_longitude(v[6])
        # vk = tools.projections.normalize_latitude(v[7])
        uk = v[6]
        vk = v[7]
        # create face's boundary
        uk_rad = np.radians(v[6])
        vk_rad = np.radians(v[7])
        XB, YB = tools.projections.boundary(U_B, V_B, R=R, uk=uk_rad, vk=vk_rad)
        point_array = arcpy.Array([arcpy.Point(x, y) for x, y in zip(XB, YB)])
        print(f'Creating face: {k}, cartographic pole, u: {v[6]}, v: {v[7]}.')
        # Combine X and Y into a list of tuples
        boundary_coords = list(zip(XB, YB))
        # Create a SpatialReference object from the WKT string
        wkt = tools.projections.update_wkt_projection(new_lon=vk, new_lat=uk)
        sr = arcpy.SpatialReference(text=wkt)
        # define map spatial reference
        map.spatialReference = sr

        # Convert geodataframe to shapefile and apply projection from wkt
        boundary_poly = arcpy.Polygon(arcpy.Array([arcpy.Point(*coords) for coords in boundary_coords]))
        centroid_geometry = boundary_poly.centroid
        boundary_ap = arcpy.management.CreateFeatureclass(out_path="in_memory", out_name="boundary", geometry_type="POLYGON", spatial_reference=sr)
        desc = arcpy.Describe(boundary_ap)
        with arcpy.da.InsertCursor(boundary_ap, ["SHAPE@"]) as cursor:
            cursor.insertRow([boundary_poly])
        boundary_layer = arcpy.management.MakeFeatureLayer(boundary_ap, "boundary_layer")[0]

        # save boundary layer to temp 
        boundary_layer_path = os.path.join("processing", f"boundary_face_{i}.shp")
        arcpy.management.CopyFeatures(boundary_layer, boundary_layer_path)
        boundary_map = arcpy.management.MakeFeatureLayer(boundary_layer_path, "boundary_map")[0]
        desc = arcpy.Describe(boundary_map)

        # MERIDIANS-PARALELS
        XM, YM, XP, YP = tools.projections.graticule(
                                                    u_min=v[2]*pi/180,
                                                    u_max=v[3]*pi/180,
                                                    v_min=v[4]*pi/180,
                                                    v_max=v[5]*pi/180,
                                                    D_u=np.radians(10),
                                                    D_v=np.radians(10),
                                                    d_u=np.radians(1),
                                                    d_v=np.radians(1), 
                                                    R=R, 
                                                    uk=uk_rad,
                                                    vk=vk_rad
                                                    )
        meridians_points = []
        parallels_points = []
        for x, y in zip(XM, YM):
            point = arcpy.Point(x, y)
            meridians_points.append(point)
        # Create an Array object from the Point objects
        meridians_array = arcpy.Array(meridians_points)

        for x,y in zip(XP, YP):
            point = arcpy.Point(x, y)
            parallels_points.append(point)
        # Create an Array object from the Point objects
        parallels_array = arcpy.Array(parallels_points)

        # Create a Polyline object from the Array objects
        meridians_polyline = arcpy.Polyline(meridians_array)
        parallels_polyline = arcpy.Polyline(parallels_array)
        arcpy.management.CopyFeatures(meridians_polyline, "in_memory/meridians_not_split")
        arcpy.management.CopyFeatures(parallels_polyline, "in_memory/parallels_not_split")
        # Split lines
        arcpy.management.SplitLine("in_memory/meridians_not_split", "in_memory/meridians_split")
        arcpy.management.SplitLine("in_memory/parallels_not_split", "in_memory/parallels_split")
        # Delete temporary not-split features
        arcpy.management.Delete("in_memory/meridians_not_split")
        arcpy.management.Delete("in_memory/parallels_not_split")
  
        # Remove lines longer than threshold using cursors
        fields = ['OID@', 'Shape_Length']
        with arcpy.da.UpdateCursor("in_memory/parallels_split", ["OID@", "SHAPE@"]) as cursor:
            for row in cursor:
                length = row[1].length
                if length > 5000000:
                    cursor.deleteRow() 

        with arcpy.da.UpdateCursor("in_memory/meridians_split", ["OID@", "SHAPE@"]) as cursor:
            for row in cursor:
                length = row[1].length
                if length > 5000000:
                    cursor.deleteRow()
    
        # create feature classes and set spatial reference
        parallels_ap = arcpy.management.CreateFeatureclass(
                out_path="in_memory",
                out_name="parallels_ap",
                geometry_type="POLYLINE",
                spatial_reference=sr
                )  
        # append features to parallels_ap
        with arcpy.da.InsertCursor(parallels_ap, ["SHAPE@"]) as cursor:
            for row in arcpy.da.SearchCursor("in_memory/parallels_split", ["SHAPE@"]):
                try:
                    cursor.insertRow(row)
                except:
                    print("Error inserting row.")
        meridians_ap = arcpy.management.CreateFeatureclass(
                out_path="in_memory",
                out_name="meridians_ap",
                geometry_type="POLYLINE",
                spatial_reference=sr
                )
        
        arcpy.management.AddField("in_memory/parallels_split", "Shape_Length", "DOUBLE")
        arcpy.management.CalculateField("in_memory/parallels_split", "Shape_Length", "!shape.length!", "PYTHON3")
        arcpy.management.AddField("in_memory/meridians_split", "Shape_Length", "DOUBLE")
        arcpy.management.CalculateField("in_memory/meridians_split", "Shape_Length", "!shape.length!", "PYTHON3")
        
        # append features to meridians_ap
        with arcpy.da.InsertCursor(meridians_ap, ["SHAPE@"]) as cursor:
            for row in arcpy.da.SearchCursor("in_memory/meridians_split", ["SHAPE@"]):
                try:
                    cursor.insertRow(row)
                except:
                    print("Error inserting row.")
        meridians_layer = arcpy.management.MakeFeatureLayer(meridians_ap, "meridians_layer")[0]
        meridians_layer_path = os.path.join("processing", f"meridians_face_no{i}.shp")
        arcpy.management.CopyFeatures(meridians_layer, meridians_layer_path)
        meridians_map = arcpy.management.MakeFeatureLayer(meridians_layer_path, "meridians_map")[0]

        parallels_layer = arcpy.management.MakeFeatureLayer(parallels_ap, "parallels_layer")[0]
        parallels_layer_path = os.path.join("processing", f"parallels_face_no{i}.shp")
        arcpy.management.CopyFeatures(parallels_layer, parallels_layer_path)
        parallels_map = arcpy.management.MakeFeatureLayer(parallels_layer_path, "parallels_map")[0]


        # MAPS       
        # 1. layers to map
        map.rotation = 60  # Nastavení rotace mapy na 60 stupňů
        map.addLayer(boundary_map)
        map.addLayer(parallels_map)
        map.addLayer(meridians_map)
        map.addDataFromPath(data_path=str(basemaps[basemap]), web_service_type="AUTOMATIC")
        boundary_layer = map.listLayers("boundary_map")[0]
        meridian_layer = map.listLayers("meridians_map")[0]
        parallels_layer = map.listLayers("parallels_map")[0]

        # set symbology
        symb_bounds = boundary_layer.symbology
        symb_bounds.renderer.symbol.color = {'RGB': [40, 40, 175, 0]}
        symb_bounds.renderer.symbol.width = 0.0
        boundary_layer.symbology = symb_bounds

        symb_meridians = meridian_layer.symbology
        symb_meridians.renderer.symbol.color = {'RGB': [0, 0, 0, 255]}
        symb_meridians.renderer.symbol.width = 0.1
        meridian_layer.symbology = symb_meridians
        parallels_layer.symbology = symb_meridians
        
        #2. pentagon
        # packing_north hemisphere
        pentagon = tools.layout.pent_create(a = side_len, angle = v[10])
        pentagon_shift = tools.layout.pent_move(pentagon, x_shift=v[8]+shift_x[i], y_shift=v[9]+shift_y[i])
        pent_center = tools.layout.pent_center(pentagon_shift)
        # print(f'{k} has center at {pent_center}')
        frame = tools.layout.frame_points(pentagon_shift)
        
        # 3. map to frame
        extent = desc.extent
        projected_extent = extent.projectAs(sr)
        map_frame = layout.createMapFrame(frame, map)
        map_frame.map = map
        # map_frame.map.addLayer(tile_layer, "BOTTOM")
        # map_frame.map.rotation = 60
        map_frame.map.spatialReference = sr
        map_frame.camera.setExtent(extent)
        map_frame.camera.heading = rotate[i]
        x_center = centroid_geometry.X
        y_center = centroid_geometry.Y
        map_frame.camera.X = x_center + face_definitions.shifts_x[i]*1000
        map_frame.camera.Y = y_center + face_definitions.shifts_y[i]*1000
        map.removeLayer(boundary_map)
    i +=1
# pro_project.activeView.refresh()
pro_project.saveACopy('test/globes_test.aprx')
layout.exportToPDF("globe_faces.pdf")

