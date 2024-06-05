import arcpy
import numpy as np
import matplotlib.pyplot as plt
import os

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
    #uv to sd
    s, d = uv_to_sd(u, v, uk, vk)

    #finding and deleting singularities
    s_min = 5 * np.pi / 180
    idx = np.where(s < s_min)[0]

    s = s[idx]
    d = d[idx]

    XC, YC = proj(R, s, d)

    return XC, YC


def graticule(umin, umax, vmin, vmax, Du, Dv, du, dv, R, uk, vk, u0, proj):
    """Creates parallels and meridians for a graticule."""
    XP = []
    YP = []
    for u in np.arange(umin, umax + Du, Du):
        vp = np.arange(vmin, vmax + dv, dv)
        up = np.full_like(vp, u)
        sp, dp = uv_to_sd(up, vp, uk, vk)
        xp, yp = proj(R, sp, dp)
        XP.append(xp)
        YP.append(yp)

    XP = np.array(XP)
    YP = np.array(YP)

    XM = []
    YM = []
    for v in np.arange(vmin, vmax + Dv, Dv):
        um = np.arange(umin, umax + du, du)
        vm = np.full_like(um, v)
        sm, dm = uv_to_sd(um, vm, uk, vk)
        xm, ym = proj(R, sm, dm)
        XM.append(xm)
        YM.append(ym)

    XM = np.array(XM)
    YM = np.array(YM)

    return XM, YM, XP, YP


def gnom(R, s, d):
    """Gnomonic projection."""
    x = R * np.tan(np.pi/2 - s) * np.cos(d)
    y = R * np.tan(np.pi/2 - s) * np.sin(d)
    return x, y


def globeface(u_min, u_max, v_min, v_max, D_u, D_v, d_u, d_v, uk, vk, R, u0, proj, ub, vb):
    plt.figure()
    plt.axis('equal')

    XM, YM, XP, YP = graticule(u_min, u_max, v_min, v_max, D_u, D_v, d_u, d_v, uk, vk, R, u0, proj)
    plt.plot(XM.T, YM.T, 'k')
    plt.plot(XP.T, YP.T, 'k')

    # Create feature class for meridians
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
    
    plt.show()

def boundary(u_min, u_max, v_min, v_max, D_u, D_v, d_u, d_v, uk, vk, R, u0, proj, ub, vb):
    plt.figure()
    plt.axis('equal')

    # Draw graticule
    XM, YM, XP, YP = graticule(u_min, u_max, v_min, v_max, D_u, D_v, d_u, d_v, uk, vk, R, u0, proj)
    plt.plot(XM.T, YM.T, 'k')
    plt.plot(XP.T, YP.T, 'k')

    # Create feature class for meridians
    meridians_output = r"E:\Work\MATLAB\meridians.shp"  
    arcpy.CreateFeatureclass_management(os.path.dirname(meridians_output), os.path.basename(meridians_output), "POLYLINE")
    
    # Add fields to the feature class
    arcpy.AddField_management("meridians.shp", "Name", "TEXT")
    
    # Create insert cursor and insert features (meridians)
    with arcpy.da.InsertCursor("meridians.shp", ["SHAPE@", "Name"]) as cursor:
        for i in range(XM.shape[0]):
            array = arcpy.Array([arcpy.Point(x, y) for x, y in zip(XM[i], YM[i])])  
            polyline = arcpy.Polyline(array)
            cursor.insertRow([polyline, f"Meridian {i + 1}"])

    # Draw continents (replace with your actual continents data and paths, using ArcPy)
    drawContinents("eur.txt", uk, vk, 6380, 0, proj)  
    # ... (Similarly for other continents and barriers)

    # Convert and project barriers
    sb, db = uv_to_sd(ub, vb, uk, vk)
    xb, yb = proj(R, sb, db, u0)

    # Draw barriers
    plt.plot(xb, yb, 'r')
    
    plt.show()
    
# ---- Main script ---- 

# Input parameters 
Du = np.radians(10)
Dv = Du
du = np.radians(1)
dv = du
R = 90
proj = gnom # Assuming you have the actual implementation of 'gnom'

# Face1 parameters
umin = np.radians(30)
umax = np.pi / 2
vmin = 0
vmax = 2 * np.pi
uk = np.pi / 2
vk = 0
u0 = 0

# Call functions
# ... (e.g., call globeface, boundary, etc. as needed)

#Example of calling function
#globeface(umin, umax, vmin, vmax, Du, Dv, du, dv, uk, vk, R, u0, proj, ub, vb) # Replace ub and vb with actual values
