"""
Skript tvoří  dvanáctistranný polyedrický glóbus složený z rovnostranných pětiúhelníků.
Parametry skriptu:
 - měřítko
 - offset od levého dolního okraje papíru
 - podkladová mapa 
 - poloměr Země
 - vrstvy pro vizualizaci (defaultně řeky + labels, kontinenty + labels)
Autor: Jakub Zapletal

Created: 2024-05-26
Last modified: 2025-06-12


"""

import arcpy
import numpy as np
import os
import re
import math
from collections import namedtuple
import time
from arcpy.cim import *

# =============================================================================
# 1. GLOBAL SETTINGS AND PARAMETERS
# =============================================================================
SCALE = 77000000
GLOBAL_OFFSET_X_MM = 50
GLOBAL_OFFSET_Y_MM = 70
PROJECT_PATH = r"project\globes.aprx"
WORKSPACE_PATH = os.path.dirname(os.path.abspath(__file__))

OUTPUT_PDF_NAME = "globe_faces.pdf"
continents_shp = r"data\_admin\continents\World_Continents.shp"
waterways_shp = (
    r"data\_physical\water\natural_earth_rivers\ne_110m_rivers_lake_centerlines.shp"
)
BASEMAP_NAME = "World_Shaded_Relief"
basemaps = {
    "Mapnik_OSM": "https://tile.openstreetmap.org/{z}/{x}/{y}.png",
    "ESRI_Imagery": "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
    "World_Imagery": "https://services.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
    "World_Street_Map": "https://services.arcgisonline.com/ArcGIS/rest/services/World_Street_Map/MapServer/tile/{z}/{y}/{x}",
    "World_Topo_Map": "https://services.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/tile/{z}/{y}/{x}",
    "World_Shaded_Relief": "https://services.arcgisonline.com/arcgis/rest/services/World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}",
    "None": None,
}
EARTH_RADIUS_M = 6378000
ua = 26.5651
ug = 52.6226
ud = 10.8123
BASE_SIDE_LENGTH = 60.0


# =============================================================================
# 2. TOOLS AND FUNCTIONS
# =============================================================================
class PentagonShifts:
    def __init__(self, d):
        self.d = d

    def shift_s(self):
        return math.sin(math.pi / 5) * self.d / math.tan(math.pi / 5)

    def shift_s2(self):
        return math.sin(2 * math.pi / 5) * self.d / math.tan(math.pi / 5)

    def shift_c(self):
        return math.cos(math.pi / 5) * self.d / math.tan(math.pi / 5)

    def shift_c2(self):
        return math.cos(2 * math.pi / 5) * self.d / math.tan(math.pi / 5)

    def shift_d(self):
        return 2 * (math.tan(0.3 * math.pi) * (self.d / 2))


def gnomonic_projection(R, s, d):
    x = R * np.tan(np.pi / 2 - s) * np.cos(d)
    y = R * np.tan(np.pi / 2 - s) * np.sin(d)
    return x, y


def uv_to_sd(u, v, uk, vk):
    dv = vk - v
    s = np.arcsin(np.sin(u) * np.sin(uk) + np.cos(u) * np.cos(uk) * np.cos(dv))
    d = -1 * np.arctan2(
        np.cos(u) * np.sin(dv),
        np.cos(u) * np.sin(uk) * np.cos(dv) - np.sin(u) * np.cos(uk),
    )
    return s, d


def graticule(u_min, u_max, v_min, v_max, D_u, D_v, d_u, d_v, R, uk, vk):
    XP, YP = [], []
    rotation_matrix = np.array([[0, 1], [-1, 0]])

    for u in np.arange(u_min, u_max + D_u, D_u):
        vp = np.arange(v_min, v_max + d_v, d_v)
        up = np.repeat(u, len(vp))
        sp, dp = uv_to_sd(up, vp, uk, vk)
        xp, yp = gnomonic_projection(R, sp, dp)
        XP.extend(xp)
        YP.extend(yp)

    XM, YM = [], []
    for v in np.arange(v_min, v_max + D_v, D_v):
        um = np.arange(u_min, u_max + d_u, d_u)
        vm = np.repeat(v, len(um))
        sm, dm = uv_to_sd(um, vm, uk, vk)
        xm, ym = gnomonic_projection(R, sm, dm)
        XM.extend(xm)
        YM.extend(ym)

    meridians = np.array([XM, YM])
    parallels = np.array([XP, YP])

    meridians_rot = np.dot(rotation_matrix, meridians)
    parallels_rot = np.dot(rotation_matrix, parallels)

    XM_rot = meridians_rot[0]
    YM_rot = meridians_rot[1]
    XP_rot = parallels_rot[0]
    YP_rot = parallels_rot[1]

    return XM_rot, YM_rot, XP_rot, YP_rot


def create_boundary_geometry(u_coords, v_coords, R, uk_rad, vk_rad):
    u_rad, v_rad = np.radians(u_coords), np.radians(v_coords)
    rotation_matrix = np.array([[0, 1], [-1, 0]])
    s_coords, d_coords = uv_to_sd(u_rad, v_rad, uk_rad, vk_rad)
    xb, yb = gnomonic_projection(R, s_coords, d_coords)
    boundary_points = np.dot(rotation_matrix, np.array([xb, yb]))
    return arcpy.Polygon(
        arcpy.Array(
            [arcpy.Point(x, y) for x, y in zip(boundary_points[0], boundary_points[1])]
        )
    )


def update_wkt_projection(new_lon, new_lat):
    wkt_string = """PROJCS["Gnomic",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984", 6378137.0, 298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Gnomonic"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Longitude_Of_Center",0.000],PARAMETER["Latitude_Of_Center",0.000],UNIT["Meter",1.0]]"""
    wkt_string = re.sub(
        r'PARAMETER\["Longitude_Of_Center",.*?\]',
        f'PARAMETER["Longitude_Of_Center",{new_lon}]',
        wkt_string,
    )
    wkt_string = re.sub(
        r'PARAMETER\["Latitude_Of_Center",.*?\]',
        f'PARAMETER["Latitude_Of_Center",{new_lat}]',
        wkt_string,
    )
    return wkt_string


def create_pentagon_frame(side_len, angle_rad, x_shift, y_shift):
    R = side_len / (2 * math.sin(math.pi / 5))
    points = [
        arcpy.Point(
            (R * math.cos(angle_rad + i * 2 * math.pi / 5)) + x_shift,
            (R * math.sin(angle_rad + i * 2 * math.pi / 5)) + y_shift,
        )
        for i in range(5)
    ]
    points.append(points[0])
    return arcpy.Polygon(arcpy.Array(points))


# =============================================================================
# 3. DEFINITION OF INDIVIDUAL FACES
# =============================================================================

FaceParams = namedtuple(
    "FaceParams",
    [
        "u_coords",
        "v_coords",
        "u_min",
        "u_max",
        "v_min",
        "v_max",
        "pole_lat",
        "pole_lon",
        "base_layout_x",
        "base_layout_y",
        "pentagon_angle_rad",
        "camera_heading",
        "camera_x_shift",
        "camera_y_shift",
    ],
)
FACES = {
    "Face 1": FaceParams(
        [-ud, ud, ug, ug, ud, -ud],
        [0, 36, 36, 324, 324, 0],
        -20,
        62,
        -41,
        41,
        (ug + ud) / 2,
        0,
        150,
        320,
        19 * math.pi / 10,
        0,
        0,
        0,
    ),
    "Face 2": FaceParams(
        [-ud, ud, ug, ug, ud, -ud],
        [72, 108, 108, 36, 36, 72],
        -20,
        62,
        31,
        112,
        (ug + ud) / 2,
        72,
        150,
        320,
        19 * math.pi / 10,
        72,
        0,
        0,
    ),
    "Face 3": FaceParams(
        [-ud, ud, ug, ug, ud, -ud],
        [144, 180, 180, 108, 108, 144],
        -20,
        62,
        105,
        185,
        (ug + ud) / 2,
        144,
        150,
        320,
        27 * math.pi / 10,
        144,
        0,
        0,
    ),
    "Face 4": FaceParams(
        [-ud, ud, ug, ug, ud, -ud],
        [216, 252, 252, 180, 180, 216],
        -20,
        62,
        175,
        260,
        (ug + ud) / 2,
        216,
        150,
        320,
        19 * math.pi / 10,
        -144,
        0,
        0,
    ),
    "Face 5": FaceParams(
        [-ud, ud, ug, ug, ud, -ud],
        [288, 324, 324, 252, 252, 288],
        -20,
        62,
        240,
        330,
        (ug + ud) / 2,
        288,
        150,
        320,
        19 * math.pi / 10,
        -72,
        0,
        0,
    ),
    "Face 6": FaceParams(
        [-ug, -ug, -ud, ud, -ud, -ug],
        [0, 72, 72, 36, 0, 0],
        -62,
        20,
        -5,
        77,
        ((-ug) + ud) / 2,
        36,
        0,
        0,
        math.pi / 10,
        0,
        0,
        0,
    ),
    "Face 7": FaceParams(
        [-ug, -ug, -ud, ud, -ud, -ug],
        [72, 144, 144, 108, 72, 72],
        -62,
        20,
        65,
        155,
        ((-ug) + ud) / 2,
        108,
        0,
        0,
        math.pi / 10,
        -72,
        0,
        0,
    ),
    "Face 8": FaceParams(
        [-ug, -ug, -ud, ud, -ud, -ug],
        [144, 216, 216, 180, 144, 144],
        -62,
        20,
        140,
        220,
        ((-ug) + ud) / 2,
        180,
        0,
        0,
        math.pi / 10,
        -144,
        0,
        0,
    ),
    "Face 9": FaceParams(
        [-ug, -ug, -ud, ud, -ud, -ug],
        [216, 288, 288, 252, 216, 216],
        -62,
        20,
        210,
        295,
        ((-ug) + ud) / 2,
        252,
        0,
        0,
        math.pi / 10,
        144,
        0,
        0,
    ),
    "Face 10": FaceParams(
        [-ug, -ug, -ud, ud, -ud, -ug],
        [288, 360, 360, 324, 288, 288],
        -62,
        20,
        -75,
        5,
        ((-ug) + ud) / 2,
        324,
        0,
        0,
        math.pi / 10,
        72,
        0,
        0,
    ),
    "Face 11": FaceParams(
        [ug, ug, ug, ug, ug, ug],
        [36, 108, 180, 252, 324, 36],
        38,
        90,
        -180,
        180,
        90,
        0,
        150,
        320,
        math.pi / 10,
        0,
        0,
        0,
    ),
    "Face 12": FaceParams(
        [-ug, -ug, -ug, -ug, -ug, -ug],
        [0, 72, 144, 216, 288, 0],
        -90,
        -38,
        -180,
        180,
        -90,
        0,
        0,
        0,
        19 * math.pi / 10,
        36,
        0,
        0,
    ),
}


# =============================================================================
# 4. MAIN
# =============================================================================


def main():
    print("Spouštím generování glóbu...")
    arcpy.env.overwriteOutput = True
    if not os.path.exists(WORKSPACE_PATH):
        print(f"Chyba: Pracovní adresář neexistuje: {WORKSPACE_PATH}")
        return
    os.chdir(WORKSPACE_PATH)
    if not os.path.exists('processing'):
        os.makedirs('processing')


    pro_project = arcpy.mp.ArcGISProject(PROJECT_PATH)
    for map_item in pro_project.listMaps():
        pro_project.deleteItem(map_item)
    for layout_item in pro_project.listLayouts():
        pro_project.deleteItem(layout_item)
    layout = pro_project.createLayout(420, 594, "MILLIMETER", "A2_Layout_Globe")
    print("Projekt a layout připraveny.")

    Shifts = PentagonShifts(BASE_SIDE_LENGTH)
    shift_x = [
        0,
        Shifts.shift_s2(),
        Shifts.shift_s(),
        -Shifts.shift_s(),
        -Shifts.shift_s2(),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ]
    shift_y = [
        -Shifts.shift_d(),
        -Shifts.shift_c2(),
        Shifts.shift_c(),
        Shifts.shift_c(),
        -Shifts.shift_c2(),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ]

    f1_params = FACES["Face 1"]
    center_x1 = f1_params.base_layout_x + shift_x[0]
    center_y1 = f1_params.base_layout_y + shift_y[0]
    shift_x[5] = center_x1 + Shifts.shift_s()
    shift_y[5] = center_y1 - Shifts.shift_c()
    shift_x[11] = center_x1 + Shifts.shift_s()
    shift_y[11] = center_y1 - Shifts.shift_c() - Shifts.shift_d()

    f12_params = FACES["Face 12"]
    center_x12 = f12_params.base_layout_x + shift_x[11]
    center_y12 = f12_params.base_layout_y + shift_y[11]
    shift_x[8] = center_x12 - Shifts.shift_s()
    shift_y[8] = center_y12 - Shifts.shift_c()
    shift_x[7] = center_x12 + Shifts.shift_s()
    shift_y[7] = center_y12 - Shifts.shift_c()
    shift_x[6] = center_x12 + Shifts.shift_s2()
    shift_y[6] = center_y12 + Shifts.shift_c2()
    shift_x[9] = center_x12 - Shifts.shift_s2()
    shift_y[9] = center_y12 + Shifts.shift_c2()

    continents_map = arcpy.management.MakeFeatureLayer(
        continents_shp, "continents_map"
    )[0]
    waterways_map = arcpy.management.MakeFeatureLayer(waterways_shp, "waterways_map")[0]

    continents_map_labels = None
    continents_shp_labels = continents_shp.replace(".shp", "_lab.shp")
    if os.path.exists(continents_shp_labels):
        continents_map_labels = arcpy.management.MakeFeatureLayer(
            continents_shp_labels, "continents_labels"
        )[0]

    symb_continents = continents_map.symbology
    if hasattr(symb_continents, "renderer"):
        symbol = symb_continents.renderer.symbol
        symbol.color = {"RGB": [128, 128, 128, 51]}
        symbol.outlineColor = {"RGB": [255, 255, 255, 255]}
        symbol.outlineWidth = 1
        continents_map.symbology = symb_continents
    if continents_map.supports("SHOWLABELS"):
        continents_map.showLabels = False
        l_cim_cont = continents_map.getDefinition("V2")
        for lc in l_cim_cont.labelClasses:
            lc.useMaplexEngine = True
            lc.expression = "$feature.CONTINENT"
            lc.maplexLabelPlacementProperties.rotationProperties.rotationExpressionInfo = arcpy.cim.CreateCIMObjectFromClassName(
                "CIMExpressionInfo", "V2"
            )
            lc.maplexLabelPlacementProperties.rotationProperties.rotationExpressionInfo.expression = (
                "90"
            )
            lc.maplexLabelPlacementProperties.rotationProperties.alignLabelToAngle = (
                False
            )
            lc.maplexLabelPlacementProperties.rotationProperties.rotationType = (
                "Arithmetic"
            )
            lc.maplexLabelPlacementProperties.rotationProperties.enable = True
            font_color = arcpy.cim.CreateCIMObjectFromClassName("CIMRGBColor", "V2")
            font_color.values = [0, 0, 0, 255]
            lc.textSymbol.symbol.color = font_color
            halo_color = arcpy.cim.CreateCIMObjectFromClassName("CIMRGBColor", "V2")
            halo_color.values = [255, 255, 255, 255]
            halo_fill = arcpy.cim.CreateCIMObjectFromClassName("CIMSolidFill", "V2")
            halo_fill.color = halo_color
            halo_symbol = arcpy.cim.CreateCIMObjectFromClassName(
                "CIMPolygonSymbol", "V2"
            )
            halo_symbol.symbolLayers = [halo_fill]

            lc.textSymbol.symbol.haloSize = 1.5
            lc.textSymbol.symbol.haloSymbol = halo_symbol
        continents_map.setDefinition(l_cim_cont)

        pro_project.save()

    symb_water = waterways_map.symbology
    if hasattr(symb_water, "renderer"):
        symb_water.renderer.symbol.color = {"RGB": [0, 0, 255, 255]}
        symb_water.renderer.symbol.width = 1
        waterways_map.symbology = symb_water
    if waterways_map.supports("SHOWLABELS"):
        waterways_map.showLabels = True
        l_cim = waterways_map.getDefinition("V2")
        for lc in l_cim.labelClasses:
            lc.expression = "$feature.name"
            font_color = arcpy.cim.CreateCIMObjectFromClassName("CIMRGBColor", "V2")
            font_color.values = (255, 105, 180, 255)
            lc.textSymbol.symbol.color = font_color
            halo_color = arcpy.cim.CreateCIMObjectFromClassName("CIMRGBColor", "V2")
            halo_color.values = [173, 216, 230, 100]
            halo_fill = arcpy.cim.CreateCIMObjectFromClassName("CIMSolidFill", "V2")
            halo_fill.color = halo_color
            halo_symbol = arcpy.cim.CreateCIMObjectFromClassName(
                "CIMPolygonSymbol", "V2"
            )
            halo_symbol.symbolLayers = [halo_fill]
            lc.textSymbol.symbol.haloSize = 1
            lc.textSymbol.symbol.haloSymbol = halo_symbol
        waterways_map.setDefinition(l_cim)

    if continents_map_labels:
        symb_labels = continents_map_labels.symbology
        if hasattr(symb_labels, "renderer"):
            symbol = symb_labels.renderer.symbol
            symbol.color = {"RGB": [255, 255, 255, 0]}
            symbol.width = 0.0
            continents_map_labels.symbology = symb_labels

        if continents_map_labels.supports("SHOWLABELS"):
            continents_map_labels.showLabels = True
            l_cim_cont = continents_map_labels.getDefinition("V2")
            for lc_cont in l_cim_cont.labelClasses:
                lc_cont.expression = "$feature.CONTINENT"
                font_color = arcpy.cim.CreateCIMObjectFromClassName("CIMRGBColor", "V2")
                font_color.values = (0, 0, 0, 255)
                lc_cont.textSymbol.symbol.color = font_color
                halo_color = arcpy.cim.CreateCIMObjectFromClassName("CIMRGBColor", "V2")
                halo_color.values = [255, 255, 255, 255]
                halo_fill = arcpy.cim.CreateCIMObjectFromClassName("CIMSolidFill", "V2")
                halo_fill.color = halo_color
                halo_symbol = arcpy.cim.CreateCIMObjectFromClassName(
                    "CIMPolygonSymbol", "V2"
                )
                halo_symbol.symbolLayers = [halo_fill]
                lc_cont.textSymbol.symbol.haloSize = 1
                lc_cont.textSymbol.symbol.haloSymbol = halo_symbol
            continents_map_labels.setDefinition(l_cim_cont)

    in_memory_fcs_to_delete = []

    for i, (face_name, params) in enumerate(FACES.items()):
        print(f"--- Zpracovávám: {face_name} ({i+1}/{len(FACES)}) ---")

        map_obj = pro_project.createMap(f"Map_{face_name.replace(' ', '_')}")
        wkt = update_wkt_projection(new_lon=params.pole_lon, new_lat=params.pole_lat)
        sr = arcpy.SpatialReference(text=wkt)
        map_obj.spatialReference = sr

        uk_rad, vk_rad = math.radians(params.pole_lat), math.radians(params.pole_lon)
        boundary_geom = create_boundary_geometry(
            params.u_coords, params.v_coords, EARTH_RADIUS_M, uk_rad, vk_rad
        )

        XM, YM, XP, YP = graticule(
            u_min=np.radians(params.u_min),
            u_max=np.radians(params.u_max),
            v_min=np.radians(params.v_min),
            v_max=np.radians(params.v_max),
            D_u=np.radians(10),
            D_v=np.radians(10),
            d_u=np.radians(1),
            d_v=np.radians(1),
            R=EARTH_RADIUS_M,
            uk=uk_rad,
            vk=vk_rad,
        )

        meridians_points = []
        parallels_points = []
        for x, y in zip(XM, YM):
            point = arcpy.Point(x, y)
            meridians_points.append(point)
        meridians_array = arcpy.Array(meridians_points)
        for x, y in zip(XP, YP):
            point = arcpy.Point(x, y)
            parallels_points.append(point)
        parallels_array = arcpy.Array(parallels_points)

        meridians_polyline = arcpy.Polyline(meridians_array)
        parallels_polyline = arcpy.Polyline(parallels_array)
        arcpy.management.CopyFeatures(
            meridians_polyline, "in_memory/meridians_not_split"
        )
        arcpy.management.CopyFeatures(
            parallels_polyline, "in_memory/parallels_not_split"
        )
        arcpy.management.SplitLine(
            "in_memory/meridians_not_split", "in_memory/meridians_split"
        )
        arcpy.management.SplitLine(
            "in_memory/parallels_not_split", "in_memory/parallels_split"
        )
        arcpy.management.Delete("in_memory/meridians_not_split")
        arcpy.management.Delete("in_memory/parallels_not_split")

        fields = ["OID@", "Shape_Length"]
        with arcpy.da.UpdateCursor(
            "in_memory/parallels_split", ["OID@", "SHAPE@"]
        ) as cursor:
            for row in cursor:
                length = row[1].length
                if length > 5000000:
                    cursor.deleteRow()

        with arcpy.da.UpdateCursor(
            "in_memory/meridians_split", ["OID@", "SHAPE@"]
        ) as cursor:
            for row in cursor:
                length = row[1].length
                if length > 5000000:
                    cursor.deleteRow()

        parallels_ap = arcpy.management.CreateFeatureclass(
            out_path="in_memory",
            out_name="parallels_ap",
            geometry_type="POLYLINE",
            spatial_reference=sr,
        )
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
            spatial_reference=sr,
        )

        arcpy.management.AddField("in_memory/parallels_split", "Shape_Length", "DOUBLE")
        arcpy.management.CalculateField(
            "in_memory/parallels_split", "Shape_Length", "!shape.length!", "PYTHON3"
        )
        arcpy.management.AddField("in_memory/meridians_split", "Shape_Length", "DOUBLE")
        arcpy.management.CalculateField(
            "in_memory/meridians_split", "Shape_Length", "!shape.length!", "PYTHON3"
        )

        with arcpy.da.InsertCursor(meridians_ap, ["SHAPE@"]) as cursor:
            for row in arcpy.da.SearchCursor("in_memory/meridians_split", ["SHAPE@"]):
                try:
                    cursor.insertRow(row)
                except:
                    print("Error inserting row.")
        meridians_layer = arcpy.management.MakeFeatureLayer(
            meridians_ap, "meridians_layer"
        )[0]
        meridians_layer_path = os.path.join("processing", f"meridians_face_no{i}.shp")
        arcpy.management.CopyFeatures(meridians_layer, meridians_layer_path)
        meridians_map = arcpy.management.MakeFeatureLayer(
            meridians_layer_path, "meridians_map"
        )[0]

        parallels_layer = arcpy.management.MakeFeatureLayer(
            parallels_ap, "parallels_layer"
        )[0]
        parallels_layer_path = os.path.join("processing", f"parallels_face_no{i}.shp")
        arcpy.management.CopyFeatures(parallels_layer, parallels_layer_path)
        parallels_map = arcpy.management.MakeFeatureLayer(
            parallels_layer_path, "parallels_map"
        )[0]

        symb_meridians = meridians_map.symbology
        if hasattr(symb_meridians, "renderer"):
            line_symbol = symb_meridians.renderer.symbol
            line_symbol.color = {"RGB": [128, 128, 128, 50]}
            line_symbol.width = 0.5
            symb_meridians.renderer.symbol = line_symbol
            meridians_map.symbology = symb_meridians

        symb_parallels = parallels_map.symbology
        if hasattr(symb_parallels, "renderer"):
            line_symbol = symb_parallels.renderer.symbol
            line_symbol.color = {"RGB": [128, 128, 128, 50]}
            line_symbol.width = 0.5
            symb_parallels.renderer.symbol = line_symbol
            parallels_map.symbology = symb_parallels

        map_obj.addLayer(waterways_map)
        map_obj.addLayer(continents_map)
        if continents_map_labels:
            map_obj.addLayer(continents_map_labels)
        map_obj.addLayer(parallels_map)
        map_obj.addLayer(meridians_map)
        for lyr in map_obj.listLayers():
            if lyr.name == "Topographic":
                map_obj.removeLayer(lyr)
                break

        if basemaps[BASEMAP_NAME] is not None:
            map_obj.addDataFromPath(
                data_path=str(basemaps[BASEMAP_NAME]), web_service_type="AUTOMATIC"
            )
        layout_x_final = params.base_layout_x + shift_x[i]
        layout_y_final = params.base_layout_y + shift_y[i]
        SCALE_FACTOR = SCALE / 77000000
        layout_x_scaled = (layout_x_final * SCALE_FACTOR) + GLOBAL_OFFSET_X_MM
        layout_y_scaled = (layout_y_final * SCALE_FACTOR) + GLOBAL_OFFSET_Y_MM
        side_length_scaled = BASE_SIDE_LENGTH * SCALE_FACTOR

        pentagon_frame_geom = create_pentagon_frame(
            side_length_scaled,
            params.pentagon_angle_rad,
            layout_x_scaled,
            layout_y_scaled,
        )
        map_frame = layout.createMapFrame(pentagon_frame_geom, map_obj)
        map_frame.name = f"MapFrame_{face_name.replace(' ', '_')}"

        extent = boundary_geom.extent
        map_frame.camera.setExtent(extent)
        map_frame.camera.heading = params.camera_heading
        map_frame.camera.X += params.camera_x_shift * SCALE_FACTOR
        map_frame.camera.Y += params.camera_y_shift * SCALE_FACTOR

    last_map_name = f"Map_{list(FACES.keys())[-1].replace(' ', '_')}"
    last_map = pro_project.listMaps(last_map_name)[0]
    for lyr in last_map.listLayers():
        print(f"{lyr.name}, IsVisible: {lyr.visible}")
    print("\nČekám 5 sekund, aby se vykreslily všechny vrstvy v mapě")
    time.sleep(5)

    output_path = os.path.join(WORKSPACE_PATH, OUTPUT_PDF_NAME)
    print(f"\nGenerování dokončeno. Exportuji do PDF: {output_path}")
    layout.exportToPDF(output_path)
    pro_project.save()

    for fc in in_memory_fcs_to_delete:
        try:
            arcpy.management.Delete(fc)
        except Exception as e:
            print(f"Chyba: {e}")

    print("Hotovo!")


if __name__ == "__main__":
    main()
