'''
Includes pixel_to_coordinates and pixel_area. pixel_to_coordinates is used inside pixel_area to fit pixel coordinates to the coordinate type from the dataset from the file. ice_olate is called from another file to provide contour and mask. The version of ice.ice_olate that will send back a mask to use is on github elligonzalez/IcebergArea.
'''

import rasterio
from rasterio.plot import show
import sys
import cv2
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import datetime
from datetime import datetime
from scipy import stats
import pandas as pd
import ice

from functools import partial    
import shapely
import shapely.ops as ops
from shapely.geometry.polygon import Polygon
from shapely.ops import orient
from shapely import wkt #geod
from pyproj import Geod, CRS
from pyproj import Proj
from pyproj import Transformer
from shapely.geometry import shape




def pixel_to_coordinates(x, y, width, height, xmin, ymin, xmax, ymax):
    ''' Function to convert pixels to coordinates.
        Args: 
                x (float): x pixel coordinate
                y (float): y pixel coordinate
                width (float): width of the image in the number of pixels, taken from the dataset
                height (float): height of the image in the number of pixels, taken from the dataset
                xmin (float): lower bound (left) of the x coordinates of the image, taken from the dataset
                ymin (float): lower bound (bottom) of the y coordinates of the image, taken from the dataset
                xmax (float): upper bound (right) of the x coordinates of the image, taken from the dataset
                ymax (float): upper bound (top) of the y coordinates of the image, taken from the dataset
        Returns:
                Converted x,y coordinates to lon, lat coordinates
    ''' 
    x_ratio = x / width
    y_ratio = y / height

    lon = xmin + (xmax - xmin) * x_ratio
    lat = ymin + (ymax - ymin) * y_ratio
    
    return lon, lat



def pixel_area(file):
    ''' Finds the area of an iceberg in km2.
        Args:
                file (string): file GEOTiff
        Return:
                Area of iceberg in km2
    '''
    areas,cont,coords,mask = ice.ice_olate(file,display=True,areaMethod="area_LL2",retCnt=True,areaThresh=7000,cutOff=8000)
    
    #Finding the interior coordinates of the contour, using mask from ice_olate    
    for contour in cont:
        coordinates = np.column_stack(np.where(mask > 0))
    
        interior_coordinates = []
        for x, y in coordinates:
            interior_coordinates.append((y,x))
    
    
    #Check interior coordinate plot
    x_coord = []
    y_coord = []
    for coordinate in interior_coordinates:
        x_coord.append(coordinate[0])
        y_coord.append(coordinate[1])
    plt.plot(x_coord, y_coord)
     
    plt.xlabel('x - axis')
    plt.ylabel('y - axis')
     
    plt.title('Interior')
    plt.show()
    

    #Finding data from image    
    dataset = rasterio.open(file)
    print(dataset.height, dataset.width)
    print(dataset.crs)
    print(dataset.bounds)
    print(dataset.transform)
    
    xmin = dataset.bounds[0]
    ymin = dataset.bounds[1]
    xmax = dataset.bounds[2]
    ymax = dataset.bounds[3]
    
    width = dataset.width
    height = dataset.height
    
    
    # Convert pixel coordinates to geographic coordinates
    coordinates = []
    for point in interior_coordinates:
        x, y = pixel_to_coordinates(point[0], point[1], width, height, xmin, ymin, xmax, ymax)
        coordinates.append((x, y))




    #Make a list of all the coordinates. In each [[all four coordinates for one center],[all four coordinates for the next center]...etc]
    dataset = rasterio.open(file)
    res = dataset.res
    corners_coordinates = []
    for coordinate in coordinates:
        individual_corner = []
        corner1 = (coordinate[0] - res[0]/2, coordinate[1] + res[1]/2)
        corner2 = (coordinate[0] + res[0]/2, coordinate[1] + res[1]/2)
        corner3 = (coordinate[0] + res[0]/2, coordinate[1] - res[1]/2)
        corner4 = (coordinate[0] - res[0]/2, coordinate[1] - res[1]/2)
        individual_corner.append(corner1)
        individual_corner.append(corner2)
        individual_corner.append(corner3)
        individual_corner.append(corner4)
        corners_coordinates.append(individual_corner)


    #Transform corners in PS to Lat/Lon
    transformer = Transformer.from_crs("epsg:3031", "epsg:4326") 
    latlon_corners = []
    for coordinate in corners_coordinates:
        individual_corner = []
        onei,onej = transformer.transform(coordinate[0][0],coordinate[0][1])
        twoi,twoj = transformer.transform(coordinate[1][0],coordinate[1][1])
        threei,threej = transformer.transform(coordinate[2][0],coordinate[2][1])
        fouri,fourj = transformer.transform(coordinate[3][0],coordinate[3][1])
        individual_corner.append([onei,onej])
        individual_corner.append([twoi,twoj])
        individual_corner.append([threei,threej])
        individual_corner.append([fouri,fourj])
        latlon_corners.append(individual_corner)

    #Finding the area in m2    
    total_area = 0
    for coordinate in latlon_corners:
        all_lat = []
        all_lon = []
        all_lat.append(coordinate[0][0])
        all_lat.append(coordinate[1][0])
        all_lat.append(coordinate[2][0])
        all_lat.append(coordinate[3][0])
        all_lon.append(coordinate[0][1])
        all_lon.append(coordinate[1][1])
        all_lon.append(coordinate[2][1])
        all_lon.append(coordinate[3][1])
        pixel = {"type": "Polygon", "coordinates": [
            [(coordinate[0][1],coordinate[0][0]),
             (coordinate[1][1],coordinate[1][0]),
             (coordinate[2][1],coordinate[2][0]),
             (coordinate[3][1],coordinate[3][0]),
             (coordinate[0][1],coordinate[0][0])]]}
        lon, lat = zip(*pixel['coordinates'][0])
        pa = Proj('+proj=aea +lat_1=' + str(min(all_lat)) + ' +lat_2=' + str(max(all_lat)) + ' +lat_0=' + str((max(all_lat)-min(all_lat))+min(all_lat)) + ' +lon_0=' + str(((max(all_lon)-min(all_lon))+min(all_lon))))
        x, y = pa(lon, lat)
        cop = {"type": "Polygon", "coordinates": [zip(x, y)]}
        total_area += shape(cop).area
    print(total_area)

    print('BERG TOTAL KM2:', total_area/1e+6)

    return total_area/1e+6
    
