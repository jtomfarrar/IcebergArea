'''
Functions same as ice_olate in ice.py, but only uses new pixel_area method. Includes thresholds, cutoffs, and also a 25% boundary from the previously recorded area
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


import pickle
    
def iceStats(dataSet):
    """Display statistics of dataset
        Args:
            dataSet (dict values): data of areas
        Return:
            None, statistics are printed out
    """
    df = pd.DataFrame(list(dataSet)) #turn dataset into list 
    Q1 = df.quantile(0.25) #find Q1
    Q3 = df.quantile(0.75) #find Q3
    med = df.median() #find median
    IQR = Q3 - Q1 #find IQR
    lower_bound = Q1 -(1.5 * IQR) #find bounds
    upper_bound = Q3 +(1.5 * IQR) 
    print('MED',med) #print
    print('LOWER',lower_bound)
    print('UPPER',upper_bound)
    print('IQR',IQR)
    print('\n')
    
def findtooLarge(dataSet,limit):
    """Find values that surpass a certain value
        Args:
            dataSet (dict): data of areas
            limit (int): value that should not be surpassed
        Return:
            List of file names that produced values that are too large
    """
    tooLarge = [i for i,j in dataSet.items() if j > limit] #list of file names w values that are above the limit
    return tooLarge #return list

def getdates(dictionary):
    """Convert strings of file names to datetime objects so they an be plotted
        Args:
            dictionary (dict): dictionary of file names w corresponding areas
        Return:
            List of sorted dates
    """
    dates = []
    chenk = []
    for i,j in dictionary.items():
        r = i[::-1]
        d = r[0:10]
        d = d[::-1]
        date_time_obj = datetime.strptime(d, '%Y-%m-%d').date()
        chenk.append((date_time_obj,j))
        chenk.sort()
    return chenk


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

    x = xmin + (xmax - xmin) * x_ratio
    y = ymin + (ymax - ymin) * y_ratio
    
    return x, y


def pixel_area(mask,directory):

    try:
                
    #Finding the interior coordinates of the contour, using mask from ice_olate    
        mask_binary = (mask > 0).astype(np.uint8)
        coordinates = np.column_stack(np.where(mask_binary))
        interior_coordinates = coordinates[:, ::-1]
    
        
    
        #Finding data from image    
        dataset = rasterio.open(directory)
        
        xmin = dataset.bounds[0]
        ymin = dataset.bounds[1]
        xmax = dataset.bounds[2]
        ymax = dataset.bounds[3]
        
        width = dataset.width
        height = dataset.height
        total_pixel = width*height
        
        # Creates an upper limit to how big the iceberg can be, or the mask of the iceberg. Edit percentage to the expected iceberg limit.
        nonZero = cv2.countNonZero(mask) #count nonzero pixels
        if nonZero >= .15*total_pixel:
                print('Too big.')
                return('Too big.')
        
        
        # Convert pixel coordinates to polar stereographic coordinates
        coordinates = []
        for point in interior_coordinates:
            x, y = pixel_to_coordinates(point[0], point[1], width, height, xmin, ymin, xmax, ymax)
            coordinates.append((x, y))
    
    
    
    
        #Make a list of all the coordinates. In each [[all four coordinates for one center],[all four coordinates for the next center]...etc]
        dataset = rasterio.open(directory)
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
        return areakm
    except:
        return('Too cloudy.')

    
    
def ice_olate(directory,layerName=None,display=False,setThresh=None,areaMethod=None,areaThresh=None,cutOff=None,retCnt = False):
    """Isolate iceberg from its surroundings. Filters out images with less than 25% white pixels and more than 75% white pixels.

        Args:
            directory (str or list): file or file directory of GEOtiffs
            layerName (str, default=None): either "red" for MODIS367 images or "blue" for MODIS721 and VIIRSM1I111 images
            display (bool, default=False): display orignial image and contours, area, and pixel area
            setThresh (int, default=None): specified threshold for thresh binary function, default value is 210
            areaMethod (str, default=None): area method that should be used can be "area_PixelImageDimensions", "area_PS", "area_PS2", "area_LL",
                "area_LL2", "area_Geographic"; pixel count is calculated if no value is provided
            areaThresh (int, default=None): the area the iceberg should be near, usually the median obtained from the first round of processing, 
                value required for reprocessing
            cutOff (int, default=None): the area the iceberg should not exceed, value required for reprocessing
            retCnt (bool, default=False): return contour of image
        Return: 
            either pixel area or true area estimate, contours and coordinates can also be returned
    """
    areas_dict = {}
    areas_list =[]
    for file in os.listdir(directory):
        name = file.rstrip('.tif')
        rimg = rasterio.open(directory+file) #open image using rasterio
        img = cv2.imread(directory+file) #open image using cv2
        if display == True: #show the original
            print(file) 
            print('ORIGINAL')
            show(rimg)
    #Isolate band
        if layerName == "red": #isolate the red channel/the first band/band 3 for MODIS367 images 
            band = rimg.read(1)         
        if layerName == "blue": #isolate the red channel/the third band/band 1 for MODIS721 images and VIIRSM1I111 images
            band = rimg.read([2,3,4])    
        elif layerName == None: #binarize true color images
            band = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        if setThresh != None: #if image is being reprocessed, the threshold can be adjusted for better results
    #1. Binary threshold
            ret,thresh4 = cv2.threshold(band,setThresh,255,cv2.THRESH_BINARY) #Eg: 200 instead of 210 makes quite a difference for some images
        else:   
            ret,thresh4 = cv2.threshold(band,210,255,cv2.THRESH_BINARY)  
    #2. Connected component analysis
        componentsNumber, labeledImage, componentStats, componentCentroids = \
        cv2.connectedComponentsWithStats(thresh4, connectivity=8)
        minArea = 800 #minimum area of blob pixels
        remainingComponentLabels = [i for i in range(1, componentsNumber) if componentStats[i][4] >= minArea] #obtain labels of blobs, 0th blob is background
        filteredImage = np.where(np.isin(labeledImage, remainingComponentLabels) == True, 255, 0).astype('uint8') #filter blobs that are too small
    #3. Dilation
        kernel = np.ones((5,5), np.uint8) 
        img_dilation = cv2.dilate(filteredImage, kernel, iterations=2) #dilate image
    #4. Blur
        blur = cv2.blur(img_dilation,(15,15))
        dataset = rasterio.open(directory+file)
        total_pixel = dataset.height*dataset.width
        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    #Parameters for skipping over cropped-out images and cloudy images.
        nonZero = cv2.countNonZero(gray) #count nonzero pixels
        if nonZero<=.25*total_pixel:
            print('Cropped out.')
            areas_dict[str(name)] =('Cropped out.')
            pass
        else:
            nonZero = cv2.countNonZero(blur) #count nonzero pixels
            if nonZero >= .9*total_pixel:
                print('Too cloudy.')
                areas_dict[str(name)] =('Too cloudy.')
                pass
            else:
            #5. Otsu Threshold
                rets, thresh = cv2.threshold(blur,0,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)
                contours, hierarchy = cv2.findContours(thresh,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_NONE) #obtain contours
            #Calculate area
                if (contours): #the image has contours, sometimes it may not
                    if display == True: #display contours
                        allContours = cv2.drawContours(img, contours, -1, (0,255,0), 3) #draw all contours on image
                        print('CONTOURED')
                        plt.imshow(cv2.cvtColor(allContours, cv2.COLOR_BGR2RGB)) #show image
                        plt.show()
                #Reprocesssing
                    if cutOff and areaThresh: #reprocess if these parameters have values
                        cnt = contours[0] #first contour
                        for i in range(len(contours)): #check for closest area among contours
                            icnt = contours[i] #current contour
                            if i == 0: #the first contour
                                h, w = img.shape[:2]
                                mask = np.zeros((h, w), np.uint8)
                                cv2.drawContours(mask, [cnt],-1, 255, -1)
                                BergArea = pixel_area(mask,directory+file)
                                if type(BergArea) == float:
                                    try:
                                        if abs(BergArea-areas_list[-1]) < .25*areas_list[-1]:
                                            howClose = abs(BergArea-areas_list[-1]) #record how close the contour area is to the area threshold
                                        else:
                                            pass
                                    except IndexError:
                                        howClose = abs(BergArea-areaThresh) #record how close the contour area is to the area threshold

                                if display == True: #display contour
                                    h, w = img.shape[:2]
                                    mask = np.zeros((h, w), np.uint8)
                                    cv2.drawContours(mask, [icnt],-1, 255, -1) #create mask
                                    res = cv2.bitwise_and(img, img, mask=mask) #apply mask to original image
                                    print('CONTOUR',i) #current contour
                                    plt.imshow(cv2.cvtColor(res, cv2.COLOR_BGR2RGB)) #show image
                                    plt.show()
                                    print('AREA',BergArea) #estimated area
                            else: #conotur is not the first one
                                h, w = img.shape[:2]
                                mask = np.zeros((h, w), np.uint8)
                                cv2.drawContours(mask, [icnt],-1, 255, -1) #create mask
                                curr_BergArea = pixel_area(mask,directory+file)
                                if type(curr_BergArea) == float:
                                    try:
                                        if abs(curr_BergArea-areas_list[-1]) < .25*areas_list[-1]:                                        
                                            if abs(curr_BergArea-areas_list[-1]) < howClose and curr_BergArea < cutOff: #if area is closer 
                                                howClose = abs(curr_BergArea-areas_list[-1]) #update how close the contour area is t
                                                BergArea = curr_BergArea #update area
                                                cnt = icnt #update contour
                                                if display == True: #show contour
                                                    h, w = img.shape[:2]
                                                    mask = np.zeros((h, w), np.uint8)
                                                    cv2.drawContours(mask, [cnt],-1, 255, -1)
                                                    res = cv2.bitwise_and(img, img, mask=mask)
                                                    print('CONTOUR',i)
                                                    plt.imshow(cv2.cvtColor(res, cv2.COLOR_BGR2RGB))
                                                    plt.show()
                                                    print('AREA',BergArea)
                                    except IndexError:
                                        if abs(curr_BergArea-areaThresh) < howClose and curr_BergArea < cutOff:
                                            howClose = abs(BergArea-areaThresh) #update howClose
                                            BergArea = curr_BergArea #update area
                                            cnt = icnt #update contour
                                            if display == True: #show contour
                                                h, w = img.shape[:2]
                                                mask = np.zeros((h, w), np.uint8)
                                                cv2.drawContours(mask, [cnt],-1, 255, -1)
                                                res = cv2.bitwise_and(img, img, mask=mask)
                                                print('CONTOUR',i)
                                                plt.imshow(cv2.cvtColor(res, cv2.COLOR_BGR2RGB))
                                                plt.show()
                                                print('AREA',BergArea)
                                   
                    #No reprocessing
                    else: 
                        cnt = max(contours, key=cv2.contourArea) #find max contour
                        h, w = img.shape[:2]
                        mask = np.zeros((h, w), np.uint8)
                        cv2.drawContours(mask, [cnt],-1, 255, -1) #create mask
                        BergArea = pixel_area(mask,directory)
                        if display == True: #display image
                            h, w = img.shape[:2]
                            mask = np.zeros((h, w), np.uint8)
                            cv2.drawContours(mask, [cnt],-1, 255, -1) #create mask
                            print('MASKED')
                            res = cv2.bitwise_and(img, img, mask=mask) #apply mask to image
                            plt.imshow(cv2.cvtColor(res, cv2.COLOR_BGR2RGB)) #show image
                            plt.show()
                            print('AREA',BergArea) #area
                try:
                    print('Final Area:',BergArea)
                    if BergArea<cutOff and abs(BergArea-areas_list[-1]) < .25*areas_list[-1]:
                        areas_dict[str(name)] = BergArea #append areas
                        areas_list.append(BergArea)
                except IndexError:
                    if BergArea<cutOff and abs(BergArea-areaThresh) < .25*areaThresh:
                        areas_dict[str(name)] = BergArea #append areas
                        areas_list.append(BergArea)
                except:
                    print('Cropped')
                    areas_dict[str(name)] ='Cropped'
    print(areas_dict)
    return areas_dict
