# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 10:42:43 2021

@author: Hector
"""
import numpy as np
import matplotlib.pyplot as plt

def xwindow(x,y,b_width):
    """
    Parameters
    ----------
    x : Spectrum x data
    y : Spectrum y data
    window : Size of x data index interval 

    Returns
    -------
    New x data, with enlarged window and corresponding normalised y data,
    which is the sum of the y data within the window i.e. 
    np.sum(y[i-window:i+window])

    """
    w = b_width # for ease of writing shorten to w
    i = w # decouple, initialising i = width for loop
    yn = list()
    xn = list()
    while i < len(y)-w: # cycle through all ys (i)        
        yi = (np.sum(y[i-w:i+w]))/(2*w) # average of window of i+-w
        xi = x[i] # reference corresponding x 
        yn.append(yi) # store averaged y
        xn.append(xi) # and corresponding x
        i += w # add w for next loop iteration

    yn = np.array(yn)
    xn = np.array(xn)  
        
    return(xn,yn)


    
    