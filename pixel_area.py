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



def pixel_area(file,layerName=None):
    ''' Finds the area of an iceberg in km2. Filters out areas less than 300km and large masks due to clouds (>10% white pixels in mask).
        Args:
                file (string): file GEOTiff
        Return:
                Area of iceberg in km2
    '''

    try:
        areas,cont,coords,mask = ice.ice_olate(file,layerName,display=True,areaMethod="area_LL2",areaThresh=5000,retCnt=True)

        #Finding the interior coordinates of the contour, using mask from ice_olate    
        mask_binary = (mask > 0).astype(np.uint8)
        coordinates = np.column_stack(np.where(mask_binary))
        interior_coordinates = coordinates[:, ::-1]

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
        total_pixel = width*height

        nonZero = cv2.countNonZero(mask) #count nonzero pixels
        if nonZero >= .1*total_pixel:
                print('Too cloudy.')
                return('Too cloudy.')
        
        
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
            max_lat = max(coordinate[0][0],coordinate[1][0],coordinate[2][0],coordinate[3][0])
            min_lat = min(coordinate[0][0],coordinate[1][0],coordinate[2][0],coordinate[3][0])
            max_lon = max(coordinate[0][1],coordinate[1][1],coordinate[2][1],coordinate[3][1])
            min_lon = min(coordinate[0][1],coordinate[1][1],coordinate[2][1],coordinate[3][1])
        
            pixel = {"type": "Polygon", "coordinates": [
                [(coordinate[0][1],coordinate[0][0]),
                 (coordinate[1][1],coordinate[1][0]),
                 (coordinate[2][1],coordinate[2][0]),
                 (coordinate[3][1],coordinate[3][0]),
                 (coordinate[0][1],coordinate[0][0])]]}
            lon, lat = zip(*pixel['coordinates'][0])
            pa = Proj('+proj=aea +lat_1=' + str(min_lat) + ' +lat_2=' + str(max_lat) + ' +lat_0=' + str((max_lat-min_lat)/2) + ' +lon_0=' + str((max_lon-min_lon)/2))
            x, y = pa(lon, lat)
            cop = {"type": "Polygon", "coordinates": [zip(x, y)]}
            total_area += shape(cop).area
        areakm = total_area/1e+6

        if areakm<=300:
            print('Too small.')
            return ('Too small.')
        else:
            print('BERG TOTAL KM2:', areakm)
            return areakm
    except:
        return('Too cloudy.')
    
