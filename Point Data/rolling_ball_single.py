# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 10:49:25 2021

@author: hlancaster
"""
import copy
import numpy as np 

def find_nearest(array, value):
    array = np.asarray(array) # assign variable 'array' as from the input
    idx = (np.abs(array - value)).argmin() # find index of the closest value
    return array[idx] # return the closest value

def spectrum_interval(x,y, lb, ub):
    lowerx = np.where(x>lb)
    upperx = np.where(x<ub)
    interx = np.intersect1d(lowerx,upperx)
    x = x[interx]  
    y = y[interx]
    return(x,y)

def remove_background(x, y, rad):

    x, y = spectrum_interval(x,y, 150, 3200)    
    x = x[::-1] # reverse order (for the rolling ball filter)
    y = y[::-1]

    yr = max(y)-min(y)
    xr = max(x)-min(x)
    y = y*(xr/yr) # make square (essential)

    # add dummy data to allow the code to catch the ends of the real data,
    # deleted later    
    dummy_xlow = np.arange(min(x)-rad,min(x)-1,5)
    dummy_ylow = np.zeros(np.shape(dummy_xlow)) + y[0] * 1.5
    x = np.concatenate((dummy_xlow, x))
    y = np.concatenate((dummy_ylow, y))
    
    dummy_xhi = np.arange(max(x)+1,max(x)+rad,5)
    dummy_yhi = np.zeros(np.shape(dummy_xhi)) + y[-1] * 1.5
    x = np.concatenate((x, dummy_xhi))
    y = np.concatenate((y, dummy_yhi))
    
    alpha = dict()
    
    for i in range(len(x)):
        
        r = rad
        
        a = x[i] # center x coord
        min_ix = np.where(x == find_nearest(x, a-r)) # (1)
        min_ix = int(min_ix[0]) # (2)
        max_ix = np.where(x == find_nearest(x, a+r))
        max_ix = int(max_ix[0])
        # (1) find nearest x that corresponds to a radius r from center
        # (2) save as integer to work with slice below
        
        # data window corresponding to radius
        xs = x[min_ix+1:max_ix] # the +1 on the lower bound skips the lowest dummy dp
        ys = y[min_ix+1:max_ix] # not having +1 on the upper bound skips the highest dummy dp
                                # this ensures the bounds don't exceed the radius
                                # stopping errors later on
        
        # increase radius, so that it overlaps window, else issues occur later
        r += 1 
    
    
        run = True
        while(run):
            unsorted = np.array((xs,ys))
            sort_done = np.array((unsorted[0][unsorted[1,:].argsort()],
                                  unsorted[1][unsorted[1,:].argsort()]))
            xsort = sort_done[0]
            ysort = sort_done[1]
            for j in range(len(ys)):
                b = ysort[j] - (r**2 - (xsort[j]-a)**2)**0.5 # center y coord
                f = b + (r**2 - (xs-a)**2)**0.5
                
                diff = ys - f
                if np.all(diff >= 0):
                    alpha[i] = np.array([xs, ys, diff, a, b], dtype=object)
                    run = False
                    break
                        
            if run == True:
                temp = np.zeros([5, len(xs)])
                temp[0,:] = temp[0,:] + xs
                temp[1,:] = temp[1,:] + np.inf
                temp[2,:] = temp[2,:] + np.inf
                temp[3,:] = temp[3,:] + np.inf
                temp[4,:] = temp[4,:] + np.inf
                alpha[i] = temp
                run = False
                break
    
        
    beta = dict() 
    beta = copy.deepcopy(alpha)  
    for k in range(len(beta)-1):
        # compares x coords
        overlap_k = np.intersect1d(beta[k][0], beta[k+1][0],
                                   assume_unique=True, return_indices=True) 
    
        diff_k = beta[k][2][overlap_k[1]]
        diff_k2 = beta[k+1][2][overlap_k[2]]
        
        beta[k+1][2][overlap_k[2]] = np.minimum(diff_k, diff_k2)
        
        
        
    t = list(beta[len(beta)-1][2]) # starting at the the RHS of the spectrum
    xtemp = list(beta[len(beta)-1][0]) # for checking, will delete
    t.reverse()
    xtemp.reverse()
    temp = list()
    for p in range(len(beta)-2, -1, -1): # going through each subsequent curve for right to left
        if beta[p][0][0] != beta[p+1][0][0]:
            t.append(beta[p][2][0]) # append the left most value 
            xtemp.append(beta[p][0][0])
        else:
            temp.append(beta[p][0][0])
    
    t = [ts + 10 for ts in t]
    
    xtemp.reverse()
    t.reverse() # reverse the spectrum
    
    xtemp = np.array(xtemp)
    t = np.array(t)
    #xtemp, t = spectrum_interval(xtemp,t, 1000, 1950)
    
    return(xtemp, t)