#!/usr/bin/env python
# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------
__author__ = 'James T. Dietrich'
__contact__ = 'james.dietrich@uni.edu'
__copyright__ = '(c) James Dietrich 2019'
__license__ = 'MIT'
__date__ = '18 June 2019'
__version__ = '1.0'
__status__ = 'initial release'
__url__ = 'https://github.com/geojames/py_streamNormalcoords'

"""
Name:           xy2sn.py
Compatibility:  Python 3.7
Description:    This program is based on an academic research article:
                
    Legleiter CJ, Kyriakidis PC. 2006. Forward and Inverse Transformations 
        between Cartesian and Channel-fitted Coordinate  Systems for Meandering 
        Rivers. Mathematical Geology 38 : 927–958. DOI: 10.1007/s11004-006-9056-6
        
                The original source code from this article was in Matlab 
                (aquired by Dietrich from Legleiter, cjl@usgs.gov). This code
                is a translation of that original Matlab code with some added
                features/tweaks.
                
                Given a digitized stream centerline and a set of input x,y
                coordinates. This software will transform the input x,y 
                coordinates into a curvalinear coordinate system based on the
                stream centerline. This stream normal coordinate system is:
                ds - a downstream distance (from the most upstream point)
                xs - a cross-stream distance (a normal (orthogonal) distace to
                                              the centerline)
                
TO RUN:
    - Fill in the blanks in the Input section
    - Run the code

DATA FORMAT:    all files are CSV inputs
                - CL_pts - X,Y points describing the centerline 
                    > Minimum columns = X, Y (caps)
                    > for spatially variable search radii extra columns can be added
                    >> variable, but symmetric - add 'R' column with distances
                    >> variable, but non-symmetric - add 'LR' and 'RR' for left and right distances
                    
                - data_pts - X,Y data points to be transformed
                    > Minimum columns = X, Y (caps)
                    > additional columns will be transfered to the output

URL:            https://github.com/geojames/py_streamNormalcoords

Requires:       scipy, numpy, pandas, matplotlib

Dev ToDo:       1) create GUI or console inputs version
                2) create function version for including in other scripts
                3) port into ArcGIS or QGIS for easy inputs

AUTHOR:         James T. Dietrich
ORGANIZATION:   University of Northern Iowa
Contact:        james.dietrich@uni.edu
Copyright:      (c) James Dietrich 2019

Permission is hereby granted, free of charge, to any person obtaining a copy of 
this software and associated documentation files (the "Software"), to deal in 
the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
of the Software, and to permit persons to whom the Software is furnished to do 
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""
#------------------------------------------------------------------------------
# Imports
import scipy as sp
import scipy.signal as signal
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import path
import matplotlib.patches as patches
import os

# ** INPUTS **

# Data In (files should be CSV)
inPath   = 'C:/path_to_your_data'
CL_pts   = 'CL_filename.csv'
data_pts = 'Data_filename.csv'

# Data Out
outPath        = 'C:/path_for_your_outputs'
out_splineCL    = 'Out_CL_filename.csv'
out_streamNorm  = 'OutData_filename.csv'

# Transformation Parameters (see readme for more info)
#   nFilt  = number of filtering iterations
#   order  = polynomial order of the filter
#   window = number of points to include in the filter window
#   nDiscr = number of segments to split the centerline into
#               - a good place to start is your total length / desired segment length
#   rMax   = maximum distance ot search for points
nFilt  = 5
order  = 3
window = 5
nDiscr = 351
rMax   = 20

# draw plots for each step (True/False)
plots = True

# end INPUTS - Run the code after this point

#--------------------------
# *** FUNCTIONS ***

# convert polar coordinates (rho,phi) to x,y coordinates
def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)
# end - pol2cart

# simple euclidian distance calculator
def euclid_dist(p1,p2):
    '''Calulates the euclidian distance between point 1 and point 2\n
        \t p1 = [x,y] array \n
        \t p2 = [x,y] array'''
    # point 1 and point 2 need to be numpy arrays first
    #  they can be multiple points (must be equal length)
    if type(p1) != np.ndarray:
        p1 = np.array([p1])
    if type(p2) != np.ndarray:
        p2 = np.array([p2])
    
    dist = np.sqrt((p2[:,0] - p1[:,0])**2 + (p2[:,1] - p1[:,1])**2)

    return dist
# end - euclid_dist
    
#--------------------------
# *** MAIN ****

# check in and out paths for trailing slash
if ('/' or "'\'") not in inPath[-1]:
    inPath = inPath + '/'
if ('/' or "'\'") not in outPath[-1]:
    outPath = outPath + '/'
if not os.path.exists(outPath):
    os.mkdir(outPath)

# Read in input data
CL_pts = pd.read_csv(inPath + CL_pts)
inCenterlinePt = np.array([CL_pts.X,CL_pts.Y])
data_pts =pd.read_csv(inPath + data_pts)

# --PART 1--
#   Smoothing the input centerline
inCenterlineFilt = inCenterlinePt

# run the Savitsky-Golay filter, iterate nFilt times
#   window size from transParam
#   polynomial order form transParam
for i in range(nFilt):
    inCenterlineFilt = signal.savgol_filter(inCenterlineFilt, 
                                               window, order)

# --PART 2-- 
#   Calculate a spline curve using the smoothed centerline
#   This allows for the line to be disctretized to the transParam value
sx = np.array(inCenterlineFilt[0,:])
sy = np.array(inCenterlineFilt[1,:])

# Calc the B-spline
#   tck = knots, B-spline coefficients, degree of the spline
#     u = values of "the parameter" (from scipy docs)
tck,u = sp.interpolate.splprep([sx,sy],s=0)

# Intrepolate the spline to a new number of points (nDiscr)
#   u_n = linspace between 0-1 (fraction of the total spline length) with th
u_n = np.linspace(0,1,nDiscr)
spline_x, spline_y = sp.interpolate.splev(u_n,tck)

# cumulative distance array, along the spline for the downstream distance 
#   calculations later
cum_dist = np.zeros([nDiscr,1])
for d in range(nDiscr-1):
    cum_dist[d+1] = cum_dist[d] + euclid_dist([spline_x[d],spline_y[d]],
                                              [spline_x[d+1],spline_y[d+1]])

# create an X,Y table of the spline centerline for saving at the end
#clOut = np.array([spline_x, spline_y])
clOut = pd.DataFrame(np.hstack((np.atleast_2d(spline_x).T, 
                                np.atleast_2d(spline_y).T, 
                                cum_dist)),columns=['X','Y','DS'])

# plot the input centerline as points and the filtered line and the spline line
if plots:    
    plt.figure()
    plt.scatter(inCenterlinePt[0,:], inCenterlinePt[1,:])
    plt.plot(inCenterlineFilt[0,:], inCenterlineFilt[1,:],'g-')
    plt.plot(spline_x, spline_y,'m--')

    n = np.arange(0,CL_pts.LR.values.size,1)
    for i, txt in enumerate(n):
        plt.annotate(txt, (inCenterlinePt[0,i], inCenterlinePt[1,i]))

# --PART 3--
# derivitives of the spline to calculate the orthogonal normal vectors
#   normals are used to construct polygons to isolate data points that are
#   near each spline segment

# calc the first derivitive of the spline
dxdt, dydt = sp.interpolate.splev(u_n,tck,der=1)

# calculate the normal vector for each nDiscr point in the spline
xNormal =  -1*(dydt/np.sqrt(dydt**2+dxdt**2))
yNormal = dxdt/np.sqrt(dydt**2+dxdt**2)

if plots:
    plt.figure()
    plt.scatter(spline_x, spline_y)
    plt.plot(spline_x, spline_y,'m-')
    plt.quiver(spline_x, spline_y, dxdt, dydt,color='g')
    plt.quiver(spline_x, spline_y, xNormal, yNormal,color='red')    
    plt.axis('equal')

# convert normal vectors to polar coords
theta = np.arctan2(yNormal, xNormal)

# extract spatially variable rMax values (if present)
# create the left bank and right bank verticies for eventual polygons
#   the verticies are rMax distance (specified in the transParams) 
#   in the normal direction

if 'R' in CL_pts.columns.values.tolist():
    print("Processing variable, but symmetric rMax")
    srMax = np.interp(np.linspace(0,CL_pts.R.values.size,nDiscr),
                      np.arange(0,CL_pts.R.values.size,1),CL_pts.R.values)
    xVertexl,yVertexl = pol2cart(srMax, theta)
    xVertexr,yVertexr = pol2cart(-1*srMax, theta)
elif 'LR' in CL_pts.columns.values.tolist():
    print("Processing variable, non-symmetric rMax")
    srMax = np.interp(np.linspace(0,CL_pts.LR.values.size,nDiscr),
                      np.arange(0,CL_pts.LR.values.size,1),CL_pts.LR.values)
    xVertexl,yVertexl = pol2cart(srMax, theta)
    
    srMax = np.interp(np.linspace(0,CL_pts.RR.values.size,nDiscr),
                      np.arange(0,CL_pts.RR.values.size,1),CL_pts.RR.values)
    xVertexr,yVertexr = pol2cart(-1*srMax, theta)
else:
    print("Processing uniform rMax")
    xVertexl,yVertexl = pol2cart(rMax, theta)
    xVertexr,yVertexr = pol2cart(-1*rMax, theta)

# add the left bank and right bank vertex values to the centerline points
#   to offest them from the centerline
poly_Lx = spline_x + xVertexl
poly_Ly = spline_y + yVertexl
poly_Rx = spline_x + xVertexr
poly_Ry = spline_y + yVertexr

if plots:
    plt.figure()
    plt.plot(spline_x, spline_y,'m--')
    plt.scatter(poly_Lx, poly_Ly, c = 'g')
    plt.scatter(poly_Rx, poly_Ry, c = 'b')

# create the polygon features (stored as matplotlib paths)
#   ploygons are left bank and right bank oriented
left_poly = []
right_poly = []

for i in range(nDiscr-1):
    l_verts = [(poly_Lx[i],poly_Ly[i]),(poly_Lx[i+1],poly_Ly[i+1]),
               (spline_x[i+1],spline_y[i+1]), (spline_x[i],spline_y[i])]
    l_patch = path.Path(l_verts)
    left_poly.append(l_patch)
    
    r_verts = [(poly_Rx[i],poly_Ry[i]),(poly_Rx[i+1],poly_Ry[i+1]),
               (spline_x[i+1],spline_y[i+1]), (spline_x[i],spline_y[i])]
    r_patch = path.Path(r_verts)
    right_poly.append(r_patch)

# plot polygons (not closed but still gives the correct visual)
if plots:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(nDiscr-1):
        l_patch = patches.PathPatch(left_poly[i],edgecolor='b', facecolor='None',zorder=1)
        r_patch = patches.PathPatch(right_poly[i],edgecolor='r', facecolor='None',zorder=1)
        ax.add_patch(l_patch)
        ax.add_patch(r_patch)
    
    ax.plot(spline_x, spline_y,'k--',lw=2,zorder=9)
    ax.scatter(data_pts.X,data_pts.Y,s=3,zorder=10)
    ax.axis('equal')

# --- PART 4 ---
# Calculate the (ds) and cross stream (xs) - the stream normal coordinates
    
# initialize stream normal coord dataframe from the input data points
xy = data_pts[['X','Y']].copy()
sn = xy.copy()
sn['ds'] = np.nan
sn['xs'] = np.nan

# evaluate the left bank SN coords
print("Evaluating Left Bank Points")

# for each polygon...
for i,poly in enumerate(left_poly):
    
    # user feedback on progress
    if i%5 == 0:
        print("   LB Segment %i of %i" %(i,len(left_poly)))
    
    # check to see what points are within the polygon, separate to new DF
    #   add the ds and xs columns (nan filled)
    inPoly = xy[poly.contains_points(xy)]
    inPoly['ds'] = np.nan
    inPoly['xs'] = np.nan
    
    # if there are points inside the polygon...
    if inPoly.shape[0] != 0:
        
        # rename the x and y columns (to avoid update confusion later)
        inPoly.rename(columns={"X": "x", "Y": "y"},inplace=True)
    
        # cross-stream (XS) and downstream (DS) coords
        #   orthogaonal distance from data point to spline line
        #   p1 = spline segment start (from poly - x = poly[3][0], y = poly[3][1])
        #   p2 = spline segment end (from poly - x = poly[2][0], y = poly[2][1])
        #   p3 = data point [X,Y]
        p1=np.array([[poly.vertices[3][0], poly.vertices[3][1]]])
        p2=np.array([[poly.vertices[2][0], poly.vertices[2][1]]])
        for j in range(inPoly.shape[0]):
            p3 = np.array([[inPoly['x'].iloc[j],inPoly['y'].iloc[j]]])
            inPoly['xs'].iloc[j] = np.linalg.norm(np.cross(p2-p1,p3-p1))/np.linalg.norm(p2-p1)
    
            # calc hypotenuse dist from spline segment start to data point
            # calc distance along spline segment inculding the cum_dist
            c_dist = euclid_dist(p1,p3)
            inPoly['ds'].iloc[j] = cum_dist[i] + np.sqrt(c_dist**2 - inPoly['xs'].iloc[j]**2)
        
        # update the sn dataframe with the new DS and XS values
        #   update keeps things organized by matching the origina dataframe index
        sn.update(inPoly)

# eval right SN coords
print("Evaluating Right Bank Points")

# for each polygon...
for i,poly in enumerate(right_poly):

    # user feedback on progress
    if i%5 == 0:
        print("   RB Segment %i of %i" %(i,len(right_poly)))
        
    # check to see what points are within the polygon, separate to new DF
    #   add the ds and xs columns (nan filled)        
    inPoly = xy[poly.contains_points(xy)]
    inPoly['ds'] = np.nan
    inPoly['xs'] = np.nan
    
    # if there are points inside the polygon...
    if inPoly.shape[0] != 0:
        
        # rename the x and y columns (to avoid update confusion later)
        inPoly.rename(columns={"X": "x", "Y": "y"},inplace=True)
    
        # cross-stream (XS) and downstream (DS) coords
        #   orthogaonal distance from data point to spline line
        #   p1 = spline segment start (from poly - x = poly[3][0], y = poly[3][1])
        #   p2 = spline segment end (from poly - x = poly[2][0], y = poly[2][1])
        #   p3 = data point [X,Y]
        p1=np.array([[poly.vertices[3][0], poly.vertices[3][1]]])
        p2=np.array([[poly.vertices[2][0], poly.vertices[2][1]]])
        for j in range(inPoly.shape[0]):
            p3 = np.array([[inPoly['x'].iloc[j],inPoly['y'].iloc[j]]])
            inPoly['xs'].iloc[j] = -1 * np.linalg.norm(np.cross(p2-p1,p3-p1))/np.linalg.norm(p2-p1)
    
            # calc hypotenuse dist from spline segment start to data point
            # calc distance along spline segment inculding the cum_dist
            c_dist = euclid_dist(p1,p3)
            inPoly['ds'].iloc[j] = cum_dist[i] + np.sqrt(c_dist**2 - inPoly['xs'].iloc[j]**2)
        
        # update the sn dataframe with the new DS and XS values
        #   update keeps things organized by matching the origina dataframe index
        sn.update(inPoly)

# --- CLEAN UP ---
        
# add attributes from input points back to the stream normal dataframe
#   drop X and Y columns because we already have them
data_pts.drop(['X','Y'],axis=1,inplace=True)
sn = pd.concat([sn,data_pts], axis=1)

# find any points that were not transformed (misisng a ds coord)
#   report to user
#   export to CSV
sn_dsNull = sn[np.isnan(sn.ds)]
sn.drop(sn[np.isnan(sn.ds)].index, inplace=True)
if sn_dsNull.shape[0] != 0:
    print("-*- WARNING - there were %i points not transformed" %(sn_dsNull.shape[0]))

# plot the final transformed data
if plots:
    plt.figure()
    plt.scatter(sn['ds'],sn['xs'],c='b',s = 5)
    #plt.scatter(val_sn['DS'],val_sn['XS'],c='r',s = 3)

# --- EXPORT ---

# export stream normal data points to outPath
sn.to_csv(outPath + out_streamNorm, index=False)

# export stream normal data points to outPath
if sn_dsNull.shape[0] != 0:
    sn_dsNull.to_csv(outPath + out_streamNorm[:-4] + '_noTrans.csv', index=False)

# save spline centerline to outPath
clOut.to_csv(outPath + out_splineCL,index=False)
