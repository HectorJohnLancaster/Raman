# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 13:02:58 2021

@author: hlancaster
"""

import numpy as np
from matplotlib import pyplot as plt, patches
import copy
import sys
import time

start_time = time.perf_counter() # start process timer

i = sample_names[0]
coord = (0,0)#(-1.28,-2.41)
xs = clean_data[i][coord][2,:]
ys = clean_data[i][coord][3,:]

# set wavenumber range
lowerx = np.where(xs>200)
upperx = np.where(xs<3200)
interx = np.intersect1d(lowerx,upperx)
xs = xs[interx][::-1]
ys = ys[interx][::-1]

yr = max(ys)-min(ys) # ys range
xr = max(xs)-min(xs) # xs range
ys = ys*(xr/yr) # make square (essential, otherwise rolling-elispoid which doesn't work)
             
alpha = dict() # initialise data dictionary

pc_count = 0 # initialise percentage counter

def find_nearest(array, value):
    array = np.asarray(array) # assign variable 'array' as from the input
    idx = (np.abs(array - value)).argmin() # find index of min value = radius limit
    return array[idx] # return the closest corresponding value 

for i in range(len(xs)):
    
    r = 499 # radius of rolling-ball, RB
    
    #============ percentage counter ==============
    pc_tot = len(xs-2*r) # denominator in % calc
    pc_count += 1 # progress counter for % calc
    pc = (pc_count*100)//(pc_tot) # % calc
    sys.stdout.write('\r' + str(pc) + '%')
    sys.stdout.flush()
    #==============================================
    
    a = xs[i] # center x coord
    min_ix = np.where(xs == find_nearest(xs, a-r)) # gets lower x index
    min_ix = int(min_ix[0])
    max_ix = np.where(xs == find_nearest(xs, a+r)) # gets upper x index
    max_ix = int(max_ix[0])
    
    x = xs[min_ix:max_ix+1] # sets x values of RB range
    y = ys[min_ix:max_ix+1] # sets y values of RB range
    
    r += 1 # makes r slightly larger than index window, 
           # required to counter rounding in 'find_nearest'

    run = True
    while(run):
        unsorted = np.array((x,y)) # initialise x,y array
                                                   # for reduced processing time
        xsort = unsorted[0][unsorted[1].argsort()] # sort x according to ascending y
        ysort = unsorted[1][unsorted[1].argsort()] # sort y according to ascending y
        for j in range(len(y)):
            # circle with center a,b and radius r, with coords x,y: r^2 = (x-a)^2 + (y-b)^2
            # a = center x coord of RB
            # for each y in the RB radius
            b = ysort[j] - (r**2 - (xsort[j]-a)**2)**0.5 # center ycoord (lower) (1)
            f = b + (r**2 - (x-a)**2)**0.5 # *(2)
            
            #(1) can be thought of as running a circle up the vertical x=a line
            #     and where it touches the spectra at x=xsort[j] that gives b
            #(2) f(x) for function of circle, giving all positive circumfrence 
            #    y values corresponding to the x values
            
            diff = y - f # subtract the circle f(x) values from the Raman y data
            if np.all(diff >= 0): # if the circle is below the data 
                # this is where the ball should be for x = a
                alpha[i] = np.array([x, y, diff, a, b], dtype='object') # save data
                run = False # move to next x
                break
                    
        if run == True: # if no solution was found for this data point
            temp = np.zeros([5, len(x)]) # initialise blank entry
            temp[0,:] = temp[0,:] + x # save relevant x data
            temp[1,:] = temp[1,:] + np.inf # set other data to np.inf
            temp[2,:] = temp[2,:] + np.inf
            temp[3,:] = temp[3,:] + np.inf
            temp[4,:] = temp[4,:] + np.inf
            alpha[i] = temp # save data
            run = False
            break

    # lowerx = np.where(x>lb)
    # upperx = np.where(x<ub)
    # interx = np.intersect1d(lowerx,upperx)
    # x = x[interx]  
    # y = y[interx]   
    
beta = dict() # initialise new dict
beta = copy.deepcopy(alpha) # make a copy to preserve original dict alpha  
for k in range(len(beta)-1): # for each solution to each x value

    lower_xrange = beta[k][0] # x values for a given RB center
    next_xrange =beta[k+1][0] # x values for the next RB center
    
    # find where these RB centers have common wavenumbers 
    overlap_k = np.intersect1d(lower_xrange, next_xrange, assume_unique=True, 
                                                          return_indices=True)

    # use overlapping indices to map corresponding solution differences
    diff_k = beta[k][2][overlap_k[1]] # for the lower RB center
    diff_k2 = beta[k+1][2][overlap_k[2]] # and the next RB center
    
    # save only the lowest values between both solutions,
    # in practice, over the full wavenumber range, this will always save
    # the data point where the circle touches the spectrum
    # like tracing the bottom of the spectrum from left to right
    beta[k+1][2][overlap_k[2]] = np.minimum(diff_k, diff_k2)
    
    #!!! what about extremeties, only roll the ball to r from each
    

# agglomerate all diff values to cover full wavenumber range with corresponding xs  
# starting at the the RHS of the spectrum, get solution difference array  
sltn = list(beta[len(beta)-1][2]) 
sltn.reverse() # it comes out as the wrong way around, so must reverse
for p in range(len(beta)-2, -1, -1): # going through each subsequent curve for right to left
    if beta[p][0][0] != beta[p+1][0][0]: 
        # as long as they don't share the same lowest data point with their predecessor 
        sltn.append(beta[p][2][0]) # append the left most diff value 

sltn.reverse() # reverse the spectrum. otherwise wrong way 

fig1 = plt.figure(figsize=(8,8))    
fig2 = plt.figure()
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)
ax1.plot(xs,sltn, label = 'r = ' + str(r))
ax1.plot(xs,ys)
ax1.legend()


# #rolling_ball check
coord = 300
xc = alpha[coord][0]
yc = alpha[coord][1]
ac = alpha[coord][3]
bc = alpha[coord][4]

y1 = bc + (r**2 - (xc-ac)**2)**0.5
y2 = bc - (r**2 - (xc-ac)**2)**0.5

ax1.plot(xc, y1, color='k') 
ax1.plot(xc, y2, color='k')

RBline = ys - sltn
ax1.plot(xs, RBline, color='k', linestyle=':')

circle1 = patches.Circle((ac,bc), radius=r, color='gold')
ax1.add_patch(circle1)

# ax1.savefig('C:/Users/hlancaster/OneDrive - University College London' + \
#                     '/PhD/Python/Data/Figures/test.pdf')

fig2 = plt.figure()
lowerx = np.where(xs>200)
upperx = np.where(xs<3200)
interx = np.intersect1d(lowerx,upperx)
ta = np.array(sltn)
ta = ta[interx]
xs2 = xs[interx]
ax2.plot(xs2, ta, linestyle='', marker='.')

plt.show()
            

# End process timer
end_time = time.perf_counter()
sys.stdout.write('\n')
print("\nScript runtime: %.2f \bs" % (end_time - start_time))
