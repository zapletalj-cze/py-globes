import arcpy
import numpy as np
import os
import sys
import pyproj
from pyproj import Proj, transform
import re
import math



class projections:
    def gnom(R, s, d):
        x = R * np.tan(np.pi/2 - s) * np.cos(d)
        y = R * np.tan(np.pi/2 - s) * np.sin(d)
        return x, y


    def uv_to_sd(u, v, uk, vk):
        # Longitude difference
        dv = vk - v  
        # Transformed latitude
        s = np.arcsin(np.sin(u) * np.sin(uk) + np.cos(u) * np.cos(uk) * np.cos(dv))
        # Transformed longitude (with quadrant adjustment using atan2)
        d = -1 * np.arctan2(np.cos(u) * np.sin(dv), np.cos(u) * np.sin(uk) * np.cos(dv) - np.sin(u) * np.cos(uk))
        return s, d


    def sd_to_uv(s, d, uk, vk):
        # Calculate longitude difference
        dv = np.arctan2(np.sin(d), np.cos(d) * np.sin(uk) + np.tan(s) * np.cos(uk))
        # Calculate transformed latitude
        u = np.arcsin(np.sin(s) * np.sin(uk) + np.cos(s) * np.cos(uk) * np.cos(d))
        # Calculate longitude
        v = vk - dv
        return u, v
        

# ---- Main script ----

# Input parameters
#Du = np.radians(10)
#Dv = Du
#du = np.radians(1)
#dv = du
#R = 90
#proj = gnom # Assuming you have the actual implementation of 'gnom'

# Face1 parameters
#umin = np.radians(30)
#umax = np.pi / 2
#vmin = 0
#vmax = 2 * np.pi
#uk = np.pi / 2
#vk = 0
#u0 = 0

# Call functions
# ... (e.g., call globeface, boundary, etc. as needed)

#Example of calling function
#globeface(umin, umax, vmin, vmax, Du, Dv, du, dv, uk, vk, R, u0, proj, ub, vb) # Replace ub and vb with actual values


# ---- Function Definitions ----

def uv_to_sd(u, v, uK, vK):
    """Converts (u, v) coordinates to (s, d) coordinates."""
    u_rad = np.radians(u)
    v_rad = np.radians(v)
    uK_rad = np.radians(uK)
    vK_rad = np.radians(vK)

    #Longitude difference
    d_v = vK_rad - v_rad
    #Catographic latitude
    s = np.arcsin(np.sin(uK_rad) * np.sin(u_rad) + np.cos(uK_rad) * np.cos(u_rad) * np.cos(d_v)) * 180 / np.pi
    #Cartographic longtitude
    num = np.cos(u_rad) * np.sin(d_v)
    denom = np.cos(u_rad) * np.sin(uK_rad) * np.cos(d_v) - np.sin(u_rad) * np.cos(uK_rad)
    d = np.arctan2(num, denom) * 180 / np.pi

    return s, d


def continents(C, R, uk, vk, u0, proj):
    """Calculates coordinates for polygon geometry in a defined projection."""
    u = C[:, 0]
    v = C[:, 1]
    #Uv to sd
    s, d = uv_to_sd(u, v, uk, vk)

    #Finding and deleting singularities
    s_min = 5 * np.pi / 180
    idx = np.where(s < s_min)[0]

    s = s[idx]
    d = d[idx]

    #Projection
    XC, YC = proj(R, s, d)

    return XC, YC


def graticule(umin, umax, vmin, vmax, Du, Dv, du, dv, R, uk, vk, u0, proj):
    """Creates parallels and meridians for a graticule."""
    XP = []
    YP = []
    for u in np.arange(umin, umax + Du, Du):
        #Create parallel
        vp = np.arange(vmin, vmax + dv, dv)
        up = np.full_like(vp, u)
        #Oblique aspect
        sp, dp = uv_to_sd(up, vp, uk, vk)
        #Project parallel
        xp, yp = proj(R, sp, dp)
        #Append row
        XP.append(xp)
        YP.append(yp)

    XP = np.array(XP)
    YP = np.array(YP)

    #Meridians
    XM = []
    YM = []
    for v in np.arange(vmin, vmax + Dv, Dv):
        #Create meridians
        um = np.arange(umin, umax + du, du)
        vm = np.full_like(um, v)
        #Oblique aspect
        sm, dm = uv_to_sd(um, vm, uk, vk)
        #Project meridian
        xm, ym = proj(R, sm, dm)
        #Append row
        XM.append(xm)
        YM.append(ym)

    XM = np.array(XM)
    YM = np.array(YM)

    return XM, YM, XP, YP


def gnom(R, u, v):
    """Gnomonic projection."""
    u_rad = (90-u)*np.pi/180
    v_rad = v*np.pi/180

    x = R * np.tan(u_rad) * np.cos(v_rad)
    y = R * np.tan(u_rad) * np.sin(v_rad)
    return x, y

    def graticule(u_min, u_max, v_min, v_max, D_u, D_v, d_u, d_v, R, uk, vk, u0, proj):
        XP, YP = [], []
        for u in np.arange(u_min, u_max + D_u, D_u):
            vp = np.arange(v_min, v_max + d_v, d_v)
            up = np.full(vp.shape, u) 
            sp, dp = projections.uv_to_sd(up, vp, uk, vk)
            xp, yp = proj(R, sp, dp)
            XP.extend(xp)
            YP.extend(yp)

        XM, YM = [], []
        for v in np.arange(v_min, v_max + D_v, D_v):
            um = np.arange(u_min, u_max + d_u, d_u)
            vm = np.full(um.shape, v)
            sm, dm = projections.uv_to_sd(um, vm, uk, vk)
            xm, ym = proj(R, sm, dm)
            XM.extend(xm)
            YM.extend(ym)

        return XM, YM, XP, YP
def globeFace(u_min, u_max, v_min, v_max, D_u, D_v, d_u, d_v, uk, vk, R, u0, proj, ub, vb):
    plt.figure()
    plt.axis('equal')

    XM, YM, XP, YP = graticule(u_min, u_max, v_min, v_max, D_u, D_v, d_u, d_v, uk, vk, R, u0, proj)
    plt.plot(XM.T, YM.T, 'k')
    plt.plot(XP.T, YP.T, 'k')

    # Create feature class for meridians
    #změnit!!!!!!
    meridians_output = r"D:\petak\Documents\škola\2022_2023\letňák\Matematická kartografie\úkol_2\meridians_shp\meridians.shp"
    arcpy.CreateFeatureclass_management(os.path.dirname(meridians_output), os.path.basename(meridians_output), "POLYLINE")

    arcpy.AddField_management("meridians.shp", "Name", "TEXT")

    with arcpy.da.InsertCursor("meridians.shp", ["SHAPE@", "Name"]) as cursor:
        for i in range(XM.shape[0]):
            array = arcpy.Array([arcpy.Point(x, y) for x, y in zip(XM[i], YM[i])])
            polyline = arcpy.Polyline(array)
            cursor.insertRow([polyline, f"Meridian {i + 1}"])


    # Draw continents and barriers (Replace with your ArcPy implementation)
    drawContinents("continents/eur.txt", uk, vk, 6380, 0, proj) 
    # ... (Similarly for other continents and barriers)
    

    def boundary(u, v, R, uk, vk):
        XB, YB = [], []
        for u_v, v_v in zip(u, v):
            s, d = projections.uv_to_sd(u_v, v_v, uk, vk)
            XB_one, YB_one = projections.gnom(R, s, d)
            XB.append(XB_one)
            YB.append(YB_one)        
        return XB, YB
    
    
    def boundary(u, v, R, uk, vk):
        XB, YB = [], []

        # Define the rotation matrix for a -90-degree rotation
        rotation_matrix = np.array([[0, 1],
                                    [-1, 0]])
        
        for u_v, v_v in zip(u, v):
            s, d = projections.uv_to_sd(u_v, v_v, uk, vk)
            XB_one, YB_one = projections.gnom(R, s, d)
            XB.append(XB_one)
            YB.append(YB_one)
        
        # Convert lists to numpy array for easier manipulation
        boundary_points = np.array([XB, YB])
        
        # Apply the rotation matrix to the boundary points
        rotated_boundary_points = np.dot(rotation_matrix, boundary_points)
        
        # Extract rotated X and Y coordinates
        rotated_X = rotated_boundary_points[0]
        rotated_Y = rotated_boundary_points[1]
        
        return rotated_X, rotated_Y
    

    def update_wkt_projection(new_lon, new_lat):
        """Updates latitude and longitude of the center in a WKT projection string."""
        wkt_string = """PROJCS["North_Pole_Gnomonic",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984", 6378137.0, 298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Gnomonic"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Longitude_Of_Center",0.000],PARAMETER["Latitude_Of_Center",0.000],UNIT["Meter",1.0]]"""
        # Update longitude 
        update_wkt = re.sub(
            r'PARAMETER\["Longitude_Of_Center",.*?\]',
            f'PARAMETER["Longitude_Of_Center",{new_lon}]', 
            wkt_string
        )
        # Update latitude
        update_wkt = re.sub(
            r'PARAMETER\["Latitude_Of_Center",.*?\]',
            f'PARAMETER["Latitude_Of_Center",{new_lat}]',
            update_wkt
        )
        return update_wkt
    
    def normalize_latitude(latitude):
        """
        Recalculates a latitude value to fit within the -90 to 90 degree range.
        """
        # Ensure latitude is within -180 to 180 range before normalization
        latitude = latitude % 180  # Normalize to 0-360 range
        while latitude > 180:
            latitude = latitude -180
        # Now, normalize to -90 to 90 range
        if -90 <= latitude <= 90:  # Already within range
            return latitude
        else:
            # Calculate the offset from the valid range
            offset = abs(latitude) - 90

            # Determine if the offset is positive or negative
            if latitude > 90:  # Northern hemisphere
                return 90 - offset
            else:  # Southern hemisphere
                return -90 + offset


    def normalize_longitude(longitude):
        """
        Recalculates a longitude value to fit within the -180 to 180 degree range.
        """
        # Check if longitude is already within range
        if -180 <= longitude <= 180:
            return longitude
        else:
            # Use modulo operator to wrap longitude around
            return ((longitude + 180) % 360) - 180


class layout:
    def pent_create(a, angle):
        """
        Based on a length (a) and angle, create a list of regular pentagram vertices. 
        The pentagram has its center at [0, 0].
        """
        R = a / (2 * math.sin(math.pi / 5))
        # Calculate the vertices of the pentagram
        points = []
        for i in range(5):
            theta = angle + i * 2 * math.pi / 5
            x = R * math.cos(theta)
            y = R * math.sin(theta)
            points.append((x, y))

        return points
    
    def pent_move(points, x_shift, y_shift):
        # move pentagram vertices by x and y
        points_moved = []
        for (x, y) in points:
            points_moved.append((x + x_shift, y + y_shift))
        return points_moved

    
    def pent_rotate(points, theta_0):
        # rotate pentagram vertices by angle
        rotated_points = []
        for (x, y) in points:
            x_rot = x * math.cos(theta_0) - y * math.sin(theta_0)
            y_rot = x * math.cos(theta_0) + y * math.sin(theta_0)
            rotated_points.append((x_rot, y_rot))
        return rotated_points
    

    
    def pent_center(vertices):
        """
        Calculate the center (centroid) of a regular pentagon given its vertices.
        """
        # Calculate center coordinates
        center_x = sum(vertex[0] for vertex in vertices) / 5
        center_y = sum(vertex[1] for vertex in vertices) / 5
        return (center_x, center_y)
    
    def frame_points(points):
        # Create an array of ArcPy Point objects from the list of points
        frame = [arcpy.Point(x, y) for x, y in points]
        # Append the first point to close the polygon loop
        frame.append(arcpy.Point(points[0][0], points[0][1]))
        # Create the polygon of the map frame boundary
        map_boundary = arcpy.Polygon(arcpy.Array(frame))
        return map_boundary


   
    
    
