#!/usr/bin/env python
# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------
__author__ = 'James T. Dietrich'
__contact__ = 'james.dietrich@austin.utexas.edu'
__copyright__ = '(c) James Dietrich 2023'
__license__ = 'MIT'
__date__ = '2023 FEB 01'
__version__ = '0.0.2'
__status__ = "major release"
__url__ = "https://github.com/geojames/NextGenRiverscapes"

"""
Name:           streamNormal.py
Compatibility:  Python 3.10
Description:    These fuctions are based on an academic research article:
                
    Legleiter CJ, Kyriakidis PC. 2006. Forward and Inverse Transformations 
    between Cartesian and Channel-fitted Coordinate  Systems for Meandering 
    Rivers. Mathematical Geology 38 : 927–958. DOI: 10.1007/s11004-006-9056-6
        
        The original source code from this article was in Matlab 
        (aquired by Dietrich from Legleiter, cjl@usgs.gov). This code
        is a translation of that original Matlab code with some added
        features/tweaks.
        
        xy2sn:
        Given a digitized stream centerline and a set of input x,y
        coordinates. This software will transform the input x,y 
        coordinates into a curvalinear coordinate system based on the
        stream centerline. 
        
        This stream normal coordinate system is:
        ds - a downstream distance (from the most upstream point)
        xs - a cross-stream distance (a normal (orthogonal) distace to
                                      the centerline)
        
        sn2xy:
        Given a set of stream normal coordinates and the centerline
        that defines the coordinate system (from xy2sn). Transforms
        ds and xs coords back to x,y coordinates in the orginal
        cartesian coordinates.
        
        NOTES:
        1. UTM or Local Metric coordinates (ie. Meters) are best.
            - Feet will work also
            - Latitude/Longitude (angular measurments) will "work", but
                will not be technically be correct.
        2. The ds and xs measurments will be in the same units as the inputs
                
TO RUN:
    - use the fuctions in this file to call the transform you want

URL:            https://github.com/geojames/py_streamNormalcoords

Requires:       scipy, numpy, pandas, matplotlib

Dev ToDo:       1) create GUI or console inputs version
                2) [CHECK] create function version for including in other scripts
                3) port into ArcGIS or QGIS for easy inputs

AUTHOR:         James T. Dietrich
ORGANIZATION:   University od Texas at Austin
Contact:        james.dietrich@austin.utexas.edu
Copyright:      (c) James Dietrich 2023

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

## helper functions

# convert polar coordinates (rho,phi) to x,y coordinates
def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)
# end - pol2cart

# simple euclidian distance calculator
def euclid_dist(p1,p2):
    '''Calulates the euclidian distance between point 1 and point 2\n
        \t p1 = [x,y] or array \n
        \t p2 = [x,y] or array'''
    # point 1 and point 2 need to be numpy arrays first
    #  they can be multiple points (must be equal length)
    if type(p1) != np.ndarray:
        p1 = np.array([p1])
    if type(p2) != np.ndarray:
        p2 = np.array([p2])
    
    dist = np.sqrt((p2[:,0] - p1[:,0])**2 + (p2[:,1] - p1[:,1])**2)

    return dist
# end - euclid_dist

## MAIN FUNTIONS
#%%
def xy2sn(x, y, cl_x, cl_y, rMax, nDiscr, verbose=True, **kwargs):
    """
    Transforms XY coordinates to a new curvilinear coord. systems based on/n
    an input stream centerline.

    Parameters
    ----------
    x : array
        X coordinates of the data points.
    y : array
        Y coordinates of the data points.
    cl_x : array
        X coordinates of the centerline
    cl_y : array
        Y coordinates of the centerline
    r_max : float or array (same length as cl_x,cl_y)
        a search radius away from the centerline to search for points
        *** SEE DOCUMENTATION ***
        > single float value - uniform search distance for all segments
        > 1-D array - allows for varing search distances downstream
        > 2-D array - allows for varing search distances downstream on both 
                        the right and left bank n X 2 (col1 = rt bank, col2 = left bank)
    nDiscr : integer
        number of segments to upsample the centerline to (0 for 1-meter spacing)
        
    verbose : boolean
        True or False (default True) - print status statements
        
    kwargs: keyword arguments = additional inputs
        to use add to the call 
          xy2sn(x, y, cl_x, cl_y, r, keyword1 = value, keyword2 = value...)
          
        nFilt : number of filtering iterations (default = 5)/n
        order : polynomial order of the filter (default = 3)/n
        window : number of points to include in the filter window (default = 5)/n
        plots : True or False (default True) - provides some simple graphs
                    during the processing /n
       
    Returns
    -------
    snOut : dataframe
        dataframe of the same length as the input XY data (columns = DS, XS)
    clOut : dataframe
        dataframe containing the smoothed and resampled centerline. Used for
        sn2xy. Save this to a CSV file!!!

    """
    
    # make up mainvars form inputs
    # centerline points
    inCenterlinePt = np.array([cl_x,cl_y])
    
    # data points
    data_pts = pd.DataFrame(np.array([x,y]).T,columns = ['X','Y'])
    
    # variables 
    nFilt  = kwargs.get("nFilt",5)
    order  = kwargs.get("order",3)
    window = kwargs.get("window",5)

    # draw plots for each step (True/False)
    plots = kwargs.get("plots",True)
    
    # check nDiscr
    #   if present use the input, if not calculate the total length of the
    #   centerline and use the rounded up total distance as nDiscr (every meter)
    #   ex. total distance = 245.6, nDiscr = 246 (~one sample per meter)
    if nDiscr == 0:
        inCL_shift = np.append(inCenterlinePt[:,1:], np.array([[np.nan],[np.nan]]),axis=1)
        d_in = euclid_dist(inCenterlinePt.T,inCL_shift.T)
        nDiscr = int(np.ceil(np.nansum(d_in))*2)
    else:
        nDiscr = int(nDiscr)
    
    # print inputs
    if verbose:
        print('INPUT VARS:')
        print("\tnFilt = %i"%nFilt)
        print("\tSpline Order = %i"%order)
        print("\tWindow Size = %i"%window)
        print("\tnDiscr = %i"%nDiscr)
        
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
        plt.title("Input Centerline Points and Smoothed CL")

        n = np.arange(0,inCenterlinePt.shape[1],1)
        for i, txt in enumerate(n):
            plt.annotate(txt, (inCenterlinePt[0,i], inCenterlinePt[1,i]))

    # --PART 3--
    # derivitives of the spline to calculate the orthogonal normal vectors
    #   normals are used to construct polygons to isolate data points that are
    #   near each spline segment

    # calc the first derivitive of the spline
    dxdt, dydt = sp.interpolate.splev(u_n,tck,der=1)

    # use the derivitive to calculate the along track angle (for use )
    clOut['phi'] = np.arctan2(dydt,dxdt)

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
        plt.title("Centerline Normal Vectors")

    # convert normal vectors to polar coords
    theta = np.arctan2(yNormal, xNormal)
    
    # extract spatially variable rMax values (if present)
    # create the left bank and right bank verticies for eventual polygons
    #   the verticies are rMax distance (specified in the transParams) 
    #   in the normal direction
    
           
    if np.size(rMax)==1:   # single, uniform rMax
        print("Processing uniform rMax")
        xVertexl,yVertexl = pol2cart(rMax, theta)
        xVertexr,yVertexr = pol2cart(-1*rMax, theta)
    elif len(rMax.shape) == 1:
        print("Processing variable, but symmetric rMax")
        srMax = np.interp(np.linspace(0,inCenterlinePt.shape[1],nDiscr),
                          np.arange(0,inCenterlinePt.shape[1],1),rMax)
        xVertexl,yVertexl = pol2cart(srMax, theta)
        xVertexr,yVertexr = pol2cart(-1*srMax, theta)
    elif len(rMax.shape) == 2:
        print("Processing variable, non-symmetric rMax")
        # Check rMAx input for columns in var. search, non-sym case
        #   transpose if rows = 2 to make columns
        if rMax.shape[0] == 2:
            rMax = rMax.T
                
        print("Processing variable, non-symmetric rMax")
        srMax = np.interp(np.linspace(0,inCenterlinePt.shape[1],nDiscr),
                          np.arange(0,inCenterlinePt.shape[1],1),rMax[:,0])
        xVertexl,yVertexl = pol2cart(srMax, theta)
        
        srMax = np.interp(np.linspace(0,inCenterlinePt.shape[1],nDiscr),
                          np.arange(0,inCenterlinePt.shape[1],1),rMax[:,1])
        xVertexr,yVertexr = pol2cart(-1*srMax, theta)
    else:
        print('ERROR - rMax has the wrong shape: ',rMax.shape)

    # add the left bank and right bank vertex values to the centerline points
    #   to offest them from the centerline
    poly_Lx = spline_x + xVertexl
    poly_Ly = spline_y + yVertexl
    poly_Rx = spline_x + xVertexr
    poly_Ry = spline_y + yVertexr

    # if plots:
    #     plt.figure()
    #     plt.plot(spline_x, spline_y,'m--')
    #     plt.scatter(poly_Lx, poly_Ly, c = 'g')
    #     plt.scatter(poly_Rx, poly_Ry, c = 'b')
        
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
        ax.set_title("Search Polygons")
        
    # --- PART 4 ---
    # Calculate the (ds) and cross stream (xs) - the stream normal coordinates
        
    # initialize stream normal coord dataframe from the input data points
    data_pts['ds'] = np.nan
    data_pts['xs'] = np.nan
    xy = data_pts[['X','Y']].copy(deep=True)

    # evaluate the left bank SN coords
    if verbose:
        print("Evaluating Left Bank Points")

    # for each polygon...
    for i,poly in enumerate(left_poly):
        
        # user feedback on progress
        if verbose:
            if i%50 == 0:
                print("   LB Segment %i of %i" %(i,len(left_poly)))
        
        # check to see what points are within the polygon, separate to new DF
        #   add the ds and xs columns (nan filled)
        inPoly = data_pts[poly.contains_points(xy)].copy(deep=True)
        inPoly['ds'] = np.nan
        inPoly['xs'] = np.nan
        
        # if there are points inside the polygon...
        if inPoly.shape[0] != 0:
        
            # cross-stream (XS) and downstream (DS) coords
            #   orthogaonal distance from data point to spline line
            #   p1 = spline segment start (from poly - x = poly[3][0], y = poly[3][1])
            #   p2 = spline segment end (from poly - x = poly[2][0], y = poly[2][1])
            #   p3 = data point [X,Y]
            p1=np.repeat(np.array([[poly.vertices[3][0], poly.vertices[3][1]]]),inPoly.shape[0],axis=0)
            p2=np.repeat(np.array([[poly.vertices[2][0], poly.vertices[2][1]]]),inPoly.shape[0],axis=0)
            
            p3 = np.array([inPoly['X'].values,inPoly['Y'].values]).T
            
            # Alt Norm Eq for vectors (np.abs(x)**2)**(1./2)
            inPoly['xs'] = (np.abs(np.cross(p2-p1,p3-p1))**2)**(1./2)/np.linalg.norm(p2-p1,axis=1)

            # calc hypotenuse dist from spline segment start to data point
            # calc distance along spline segment inculding the cum_dist
            c_dist = euclid_dist(p1,p3)
            inPoly['ds'] = cum_dist[i] + np.sqrt(c_dist**2 - inPoly['xs']**2)
            
            # update the datapoints dataframe with the new DS and XS values
            #   update keeps things organized by matching the original dataframe index
            data_pts.update(inPoly)

    # eval right SN coords
    if verbose:
        print("Evaluating Right Bank Points")

    # for each polygon...
    for i,poly in enumerate(right_poly):

        # user feedback on progress
        if verbose:
            if i%50 == 0:
                print("   RB Segment %i of %i" %(i,len(right_poly)))
            
        # check to see what points are within the polygon, separate to new DF
        #   add the ds and xs columns (nan filled)        
        inPoly = data_pts[poly.contains_points(xy)].copy(deep=True)
        inPoly['ds'] = np.nan
        inPoly['xs'] = np.nan
        
        # if there are points inside the polygon...
        if inPoly.shape[0] != 0:
            
            # cross-stream (XS) and downstream (DS) coords
            #   orthogaonal distance from data point to spline line
            #   p1 = spline segment start (from poly - x = poly[3][0], y = poly[3][1])
            #   p2 = spline segment end (from poly - x = poly[2][0], y = poly[2][1])
            #   p3 = data point [X,Y]
            p1=np.repeat(np.array([[poly.vertices[3][0], poly.vertices[3][1]]]),inPoly.shape[0],axis=0)
            p2=np.repeat(np.array([[poly.vertices[2][0], poly.vertices[2][1]]]),inPoly.shape[0],axis=0)
            
            p3 = np.array([inPoly['X'].values,inPoly['Y'].values]).T
            
            # Alt Norm Eq for vectors (np.abs(x)**2)**(1./2)
            inPoly['xs'] = -1 * (np.abs(np.cross(p2-p1,p3-p1))**2)**(1./2)/np.linalg.norm(p2-p1,axis=1)

            # calc hypotenuse dist from spline segment start to data point
            # calc distance along spline segment inculding the cum_dist
            c_dist = euclid_dist(p1,p3)
            inPoly['ds'] = cum_dist[i] + np.sqrt(c_dist**2 - inPoly['xs']**2)
            
            # update the sn dataframe with the new DS and XS values
            #   update keeps things organized by matching the origina dataframe index
            data_pts.update(inPoly)

    # --- CLEAN UP ---

    # find any points that were not transformed (misisng a ds coord)
    #   report to user
    #   export to CSV
    data_dsNull = data_pts[np.isnan(data_pts.ds)]
    if data_dsNull.shape[0] != 0:
        if verbose:
            print("-*- WARNING - there were %i points not transformed -*-" %(data_dsNull.shape[0]))

    # plot the final transformed data
    if plots:
        plt.figure()
        plt.scatter(data_pts['ds'],data_pts['xs'],c='b',s = 5)
        plt.xlabel("Downstream Distance (ds)")
        plt.ylabel("Crossstream Distance (xs)")
        plt.title("Stream Normal Points")
        
    return data_pts, clOut


#%%
def sn2xy(ds, xs, cl, verbose=True, finite=False):
    """ 
    translates data points that are in a known stream-normal coordinate 
    system (eg. a known centerline processed with xy2sn.py) back into 
    cartesian XY coordinates to match the original data you processed with xy2sn.
    
    Parameters
    ----------
    ds : array (n X 1)
        Downstream coordinates of the input data.
    xs : array (n X 1)
        Crossstream coordinates of the input data.
    cl : dataframe
        Centerline output form xy2sn (Columns = X, Y, DS, phi)
    finite : boolean
        True - returns Finite Difference coorindates (from Leg. & Kyr. implementation),
        False (default)
    

    Returns
    -------
    xy_out : dataframe with transformed coordinates, same order as inputs

    """

    # ** INPUTS **
    # make a DF from downstream and crossstream coordinates
    
    sn = pd.DataFrame(np.array([ds,xs]).T,
                      columns=['ds','xs'])
    
    # extract data point with null ds/xs coords
    sn_nan = sn[sn.isna().any(axis=1)]
    sn.drop(sn[np.isnan(sn.ds)].index, inplace=True)
    
    # - add columns for X and Y transofrmed coords
    # - if finite distance option - add finite columns
    sn['xout'] = np.nan 
    sn['yout'] = np.nan 
    if finite == True:
        sn['xout_fd'] = np.nan 
        sn['yout_fd'] = np.nan 
        
    
    if verbose:
        print('SN: %i  | SN_nan: %i'%(sn.shape[0],sn_nan.shape[0]))

    for i,row in sn.iterrows():
        
        # Calculate the downstream distance difference between the target 
        #   point (row) and the centerline points. (Signed and Absolute)
        # Sort by the centeline DF by downstream distance to the target point (row)    
        
        cl_sort = cl.copy(deep=True)
        cl_sort['ds_diff'] = cl_sort['DS'] - row['ds']
        cl_sort['ds_diffABS'] = np.abs(cl_sort['DS'] - row['ds'])
        cl_sort.sort_values(by='ds_diffABS',inplace=True)
        
        # Centerline Angle Method
        # This method used the original centerline angles calculated by xy2sn
        #   conversion
        # Fine the closest 2 upstream points (top two rows of the sorted DF)
        # Sort those based on the signed diffference (neg = upstream, pos = downstream)
        # get the alongstream angle value from the closest upstream point (phi)
        cl_top = cl_sort.iloc[0:2].sort_values(by='ds_diff')
        cl_pt_phi = cl_top.phi.iloc[0]
        
        # get the abs. downstream distance difference between the target and CL point
        #   calculate a point rotation (sn,xs) > (x,y) around a 
        #   central axis (CL point)
        #       x' = x cos(phi) - y sin(phi)
        #       y' = x sin(phi) + y cos(phi)
        #   x',y' = cartesian offset coordinates
        #   x = downstream distance difference, y = cross-stream (xs) distance 
        dds = cl_top.ds_diffABS.iloc[0]
        dx = dds * np.cos(cl_pt_phi) - row['xs'] * np.sin(cl_pt_phi)
        dy = dds * np.sin(cl_pt_phi) + row['xs'] * np.cos(cl_pt_phi)
        
        # add the rotated target point coordinates to the reference CL point to
        #   get the final output x,y coordinates
        sn['xout'][i] = cl_top.X.iloc[0] + dx
        sn['yout'][i] = cl_top.Y.iloc[0] + dy
        
        # FINITE DIFFERECE
        # ** Original Method from Legleiter and Kyriakidis **
        # We want the segment bracketed by the two nodes closest to the point to 
        #   be transformed
        # Compute finite differences dx0/ds and dy0/ds    
        # The basic idea is to calculate dx0/ds and dy0/ds using a finite
        # difference approach and then to calculate the x,y coordinates of a
        # pair of s,n coordinates using the equations given in the caption of
        # Figure 1 of Smith and McLean (1984)
        if finite == True:
            dxy = cl_sort.iloc[0] - cl_sort.iloc[1]
            dy_ds = dxy.Y/dxy.ds_diff
            dx_ds = dxy.X/dxy.ds_diff
            sn['xout_fd'][i] = cl_sort.iloc[0]['X'] - row['xs'] * dy_ds
            sn['yout_fd'][i] = cl_sort.iloc[0]['Y'] + row['xs'] * dx_ds
        
    xy_out = pd.concat([sn,sn_nan])
    xy_out.sort_index(inplace=True)
        
    return xy_out