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

def area_PixelImageDimensions(mask,img):
    """Calculate area of images in Antarctic projection from pixel count and image dimensions
        Args:
            mask (array): array of pixel values obtained from image processing
            img (array): original array before processing 
        Return:
            Iceberg area in km
    """
    height, width, channels = img.shape #image pixel height and width
    area = cv2.countNonZero(mask)*(250*250) #count number of pixels in processed images, multiply by image dimensions (300x300(km)), and divide by pixel width and height
    return area #return area

def area_PS(cont,file,retCon = False):
    """Calculate area of images in Antarctic projection using polar stereographic coordinates
        Args:
            cont (array): pixel coordinates of contour
            file (string): name of file 
            retCon (bool, default=False): to return transformed coordinates
        Return:
            Iceberg area in km
    """
    with rasterio.open(file, 'r') as r: #open file 
        T0 = r.transform #transformation matrix
    temp = [] #list of coordinates
    for i in cont: #transform every contour coordinate to polar streographic coordinates
        temp.append(T0*(i[0][1],i[0][0])) #use transformation matrix for conversion
    geom = Polygon(temp) #create polygon out of converted coordinates
    if retCon == True: #return the area and the transformed coordinates
        return geom.area/1e+6,temp
    else: #return the area
        return geom.area/1e+6

def area_PS2(cont,file,retCon = False):
    """Calculate area of images in Antarctic projection using polar stereographic coordinates
        Args:
            cont (array): pixel coordinates of contour
            file (string): name of file 
            retCon (bool, default=False): to return transformed coordinates
        Return:
            Iceberg area in km
        Area estimates for this function and area_PS are the same besides the last two digits
    """
    with rasterio.open(file, 'r') as r: #open file 
        T0 = r.transform #transformation matrix
    temp = [] #list of coordinates
    for i in cont: #transform every contour coordinate to polar streographic coordinates
        temp.append(T0*(i[0][0],i[0][1])) #use transformation matrix for conversion
    geom = Polygon(temp) #create polygon out of converted coordinates
    if retCon == True: #return the area and the transformed coordinates
        return geom.area/1e+6,temp
    else: #return the area
        return geom.area/1e+6
def area_LL(cont,file,retCon =  False):
    """Calculate area of images in Antarctic projection using lat lon coordinates
        Args:
            cont (array): pixel coordinates of contour
            file (string): name of file 
            retCon (bool, default=False): to return transformed coordinates
        Return:
            Iceberg area in km
    """
    with rasterio.open(file, 'r') as r: #open file
        T0 = r.transform #transformation matrix
    temp = [] #list of polar stereographic coordinates
    for i in cont: #transform every contour coordinate to polar streographic coordinates
        temp.append(T0*(i[0][1],i[0][0])) #use transformation matrix for conversion
    ps_to_latlon = [] #list of latlon coordinates
    transformer = Transformer.from_crs("epsg:3031","epsg:4326") #epsg code for conversion
    for i in temp: #transform every ps coordinate to latlon
        (x,y) = transformer.transform(i[0],i[1])
        ps_to_latlon.append((x,y)) 
    geom = Polygon(ps_to_latlon) #create polygon out of converted coordinates
    geod = Geod(ellps="WGS84") #use geod to calculate geodesic area
    poly_area, poly_perimeter = geod.geometry_area_perimeter(orient(geom)) #caculate area and perimeter
    if retCon == True: #return the area and the transformed coordinates
        return poly_area/1e+6,ps_to_latlon
    else: #return the area
        return poly_area/1e+6
    
def area_LL2(cont,file,retCon =  False):
    """Calculate area of images in Antarctic projection using flipped lat lon coordinates
        Args:
            cont (array): pixel coordinates of contour
            file (string): name of file 
            retCon (bool, default=False): to return transformed coordinates
        Return:
            Iceberg area in km
    """
    with rasterio.open(file, 'r') as r: #open file
        T0 = r.transform #transformation matrix
    temp = [] #list of polar stereographic coordinates
    for i in cont: #transform every contour coordinate to polar streographic coordinates
        temp.append(T0*(i[0][0],i[0][1])) #use transformation matrix for conversion
    ps_to_latlon = [] #list of latlon coordinates
    transformer = Transformer.from_crs("epsg:3031","epsg:4326") #epsg code for conversion
    for i in temp: #transform every ps coordinate to latlon
        (x,y) = transformer.transform(i[0],i[1])
        ps_to_latlon.append((x,y)) 
    geom = Polygon(ps_to_latlon) #create polygon out of converted coordinates
    geod = Geod(ellps="WGS84") #use geod to calculate geodesic area
    poly_area, poly_perimeter = geod.geometry_area_perimeter(orient(geom)) #caculate area and perimeter
    if retCon == True: #return the area and the transformed coordinates
        return poly_area/1e+6,ps_to_latlon
    else: #return the area
        return poly_area/1e+6
    
def area_Geographic(cont,file,retCon =  False): #geographic projection images
    """Calculate area of geographically projected images using lat lon coordinates
        Args:
            cont (array): pixel coordinates of contour
            file (string): name of file 
            retCon (bool, default=False): to return transformed coordinates
        Return:
            Iceberg area in km
    """
    with rasterio.open(file, 'r') as r: #open file
        T0 = r.transform #transformation matrix
    temp = [] #list of latlon coordinates
    for i in cont: #transform every contour coordinate to latlon coordinates
        temp.append(T0*(i[0][1],i[0][0])) #use transformation matrix for conversion
    geom = Polygon(temp) #create polygon out of converted coordinates
    geod = Geod(ellps="WGS84") #use geod to calculate geodesic area
    poly_area, poly_perimeter = geod.geometry_area_perimeter(orient(geom)) #caculate area and perimeter
    if retCon == True: #return the area and the transformed coordinates
        return poly_area/1e+6,temp
    else: #return the area
        return poly_area/1e+6
    
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
    if areaThresh and cutOff and type(directory) != str: #during reprocessing, directory is a dictionary instead of a string
        adir = directory
    else: #directory is a string 
        adir = glob.glob(directory) 
    trackAreas = {} #dictionary of file names and areas   
    for file in adir: #loop through files in directory
        rimg = rasterio.open(file) #open image using rasterio
        img = cv2.imread(file) #open image using cv2
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
        dataset = rasterio.open(file)
        total_pixel = dataset.height*dataset.width
        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    #Parameters for skipping over cropped-out images and cloudy images.
        nonZero = cv2.countNonZero(gray) #count nonzero pixels
        if nonZero<=.25*total_pixel:
            print('Cropped out.')
            return('Cropped out.')
        else:
            nonZero = cv2.countNonZero(blur) #count nonzero pixels
            if nonZero >= .9*total_pixel:
                print('Too cloudy.')
                return('Too cloudy.')
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
                                if areaMethod=="area_PixelImageDimensions": #use first area method
                                    if retCnt==True: #error if trying to obtain contour coordinates 
                                        raise Exception('Cannot obtain contour coordinates with this area method. Please try a different area method.')
                                    h, w = img.shape[:2]
                                    mask = np.zeros((h, w), np.uint8)
                                    cv2.drawContours(mask, [cnt],-1, 255, -1) #create mask
                                    BergArea = area_PixelImageDimensions(mask,img)
                                elif areaMethod=="area_PS": #use second area method
                                    if retCnt == True: #return coordinates of contour and area
                                        BergArea,coords = area_PS(cnt,file,retCon=True)
                                    else: #return only area
                                        BergArea = area_PS(cnt,file)
                                elif areaMethod=="area_PS2": #use third area method
                                    if retCnt == True: 
                                        BergArea,coords = area_PS2(cnt,file,retCon=True)
                                    else:
                                        BergArea = area_PS2(cnt,file)
                                elif areaMethod=="area_LL": #use fourth area method
                                    if retCnt == True:
                                        BergArea,coords = area_LL(cnt,file,retCon=True)
                                    else:
                                        BergArea = area_LL(cnt,file)
                                elif areaMethod=="area_LL2": #use fifth area method
                                    if retCnt == True:
                                        BergArea,coords = area_LL2(cnt,file,retCon=True)
                                    else:
                                        BergArea = area_LL2(cnt,file)
                                elif areaMethod=="area_Geographic": #use sixth area method
                                    if retCnt == True:
                                        BergArea,coords = area_Geographic(cnt,file,retCon=True)
                                    else:
                                        BergArea = area_Geographic(cnt,file)                   
                                else: #use pixel area method
                                    if retCnt==True: #error if trying to obtain contour coordinates 
                                        raise Exception('Cannot obtain contour coordinates with this area method. Please try a different area method.')
                                    h, w = img.shape[:2]
                                    mask = np.zeros((h, w), np.uint8)
                                    cv2.drawContours(mask, [cnt],-1, 255, -1) #create mask
                                    BergArea = cv2.countNonZero(mask) #count nonzero pixels
                                howClose = abs(BergArea-areaThresh) #record how close the contour area is to the area threshold
                                if display == True: #display contour
                                    h, w = img.shape[:2]
                                    mask = np.zeros((h, w), np.uint8)
                                    cv2.drawContours(mask, [icnt],-1, 255, -1) #create mask
                                    res = cv2.bitwise_and(img, img, mask=mask) #apply mask to original image
                                    print('CONTOUR',i) #current contour
                                    plt.imshow(cv2.cvtColor(res, cv2.COLOR_BGR2RGB)) #show image
                                    plt.show()
                                    #print('AREA',BergArea) #estimated area
                            else: #conotur is not the first one
                                if areaMethod=="area_PixelImageDimensions": #use first area method
                                    if retCnt==True: #error if trying to obtain contour coordinates 
                                        raise Exception('Cannot obtain contour coordinates with this area method. Please try a different area method.')
                                    h, w = img.shape[:2]
                                    mask = np.zeros((h, w), np.uint8)
                                    cv2.drawContours(mask, [icnt],-1, 255, -1) #create mask
                                    curr_BergArea = area_PixelImageDimensions(mask,img)
                                elif areaMethod=="area_PS": #use second area method
                                    if retCnt == True: #return coordinates of contour and area
                                        curr_BergArea,curr_coords = area_PS(icnt,file,retCon=True)
                                    else: #return only area
                                        curr_BergArea = area_PS(icnt,file)
                                elif areaMethod=="area_PS2": #use third area method
                                    if retCnt == True: 
                                        curr_BergArea,curr_coords = area_PS2(icnt,file,retCon=True)
                                    else:
                                        curr_BergArea = area_PS2(icnt,file)
                                elif areaMethod=="area_LL": #use fourth area method
                                    if retCnt == True:
                                        curr_BergArea,curr_coords = area_LL(icnt,file,retCon=True)
                                    else:
                                        curr_BergArea = area_LL(icnt,file)
                                elif areaMethod=="area_LL2": #use fifth area method
                                    if retCnt == True:
                                        curr_BergArea,curr_coords = area_LL2(icnt,file,retCon=True)
                                    else:
                                        curr_BergArea = area_LL2(icnt,file)
                                elif areaMethod=="area_Geographic": #use sixth area method
                                    if retCnt == True:
                                        curr_BergArea,curr_coords = area_Geographic(icnt,file,retCon=True)
                                    else:
                                        curr_BergArea = area_Geographic(icnt,file)                   
                                else: #use pixel area method
                                    if retCnt==True: #error if trying to obtain contour coordinates 
                                        raise Exception('Cannot obtain contour coordinates with this area method. Please try a different area method.')
                                    h, w = img.shape[:2]
                                    mask = np.zeros((h, w), np.uint8)
                                    cv2.drawContours(mask, [icnt],-1, 255, -1) #create mask
                                    curr_BergArea = cv2.countNonZero(mask) #count nonzero pixels  
                                if abs(curr_BergArea-areaThresh) < howClose and curr_BergArea < cutOff: #if area is closer to estimate
                                    BergArea = curr_BergArea #update area
                                    if retCnt == True: #if coordinates are needed then update coordinates
                                        coords=curr_coords
                                    howClose = abs(BergArea-areaThresh) #update how close the contour area is to the area threshold
                                    cnt = icnt #update contour
                                    if display == True: #show contour
                                        h, w = img.shape[:2]
                                        mask = np.zeros((h, w), np.uint8)
                                        cv2.drawContours(mask, [cnt],-1, 255, -1)
                                        res = cv2.bitwise_and(img, img, mask=mask)
                                        print('CONTOUR',i)
                                        plt.imshow(cv2.cvtColor(res, cv2.COLOR_BGR2RGB))
                                        plt.show()
                                        #print('AREA',BergArea)
                    #No reprocessing
                    else: 
                        cnt = max(contours, key=cv2.contourArea) #find max contour
                        if areaMethod=="area_PixelImageDimensions": #use first area method
                            if retCnt==True: #error if trying to obtain contour coordinates 
                                raise Exception('Cannot obtain contour coordinates with this area method. Please try a different area method.')
                            h, w = img.shape[:2]
                            mask = np.zeros((h, w), np.uint8)
                            cv2.drawContours(mask, [cnt],-1, 255, -1) #create mask
                            BergArea = area_PixelImageDimensions(mask,img)
                        elif areaMethod=="area_PS": #use second area method
                            if retCnt == True: #return coordinates of contour and area
                                BergArea,coords = area_PS(cnt,file,retCon=True)
                            else: #return only area
                                BergArea = area_PS(cnt,file)
                        elif areaMethod=="area_PS2": #use third area method
                            if retCnt == True: 
                                BergArea,coords = area_PS2(cnt,file,retCon=True)
                            else:
                                BergArea = area_PS2(cnt,file)
                        elif areaMethod=="area_LL": #use fourth area method
                            if retCnt == True:
                                BergArea,coords = area_LL(cnt,file,retCon=True)
                            else:
                                BergArea = area_LL(cnt,file)
                        elif areaMethod=="area_LL2": #use fifth area method
                            if retCnt == True:
                                BergArea,coords = area_LL2(cnt,file,retCon=True)
                            else:
                                BergArea = area_LL2(cnt,file)
                        elif areaMethod=="area_Geographic": #use sixth area method
                            if retCnt == True:
                                BergArea,coords = area_Geographic(cnt,file,retCon=True)
                            else:
                                BergArea = area_Geographic(cnt,file)                   
                        else: #use pixel area method
                            if retCnt==True: #error if trying to obtain contour coordinates 
                                raise Exception('Cannot obtain contour coordinates with this area method. Please try a different area method.')
                            h, w = img.shape[:2]
                            mask = np.zeros((h, w), np.uint8)
                            cv2.drawContours(mask, [cnt],-1, 255, -1) #create mask
                            BergArea = cv2.countNonZero(mask) #count nonzero pixels
                            print(BergArea)
                        if display == True: #display image
                            h, w = img.shape[:2]
                            mask = np.zeros((h, w), np.uint8)
                            cv2.drawContours(mask, [cnt],-1, 255, -1) #create mask
                            print('MASKED')
                            res = cv2.bitwise_and(img, img, mask=mask) #apply mask to image
                            plt.imshow(cv2.cvtColor(res, cv2.COLOR_BGR2RGB)) #show image
                            plt.show()
                            #print('AREA',BergArea) #area
                    trackAreas[file] = BergArea #record area  
            #Return
            if retCnt == True:
                return trackAreas,cnt,coords,mask #return areas, pixel contours, and coordinates
            else:
                return trackAreas #return areas
