# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 10:49:25 2021

@author: hlancaster
"""
import copy
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import sys
import time
import xwindow
np.seterr(all = "raise")
start_time = time.perf_counter() 

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

import winsound
frequency = 500  # Set Frequency To 2500 Hertz
duration = 1000  # Set Duration To 1000 ms == 1 second
   
# xs = xs[~np.isnan(xs)] # remove any NaNs 
# ys = ys[~np.isnan(ys)]
# xsafe = clean_data['graphite 785'][0,:]
# ysafe = clean_data['graphite 785'][1,:]
# xsafe = sliced_data[sample_names[1]][(-100, -50)][:,2]
# ysafe = sliced_data[sample_names[1]][(-100, -50)][:,3]
# xsafe = raw_data[sample_names[0]][:,0]
# ysafe = raw_data[sample_names[0]][:,1]
# xsafe = clean_data[sample_names[-1]][0,:]
# ysafe = clean_data[sample_names[-1]][1,:]
    
rad = 3000#200
# xs = xn
# ys = yn  
xs = xsafe
ys = ysafe
lb = int(min(xs)-1) # + 200
ub = int(max(xs)+1) # - 400

# xn, yn = xwindow.xwindow(xs, ys, 2)
# xs = xn
# ys = yn

xs, ys = spectrum_interval(xs, ys, lb, ub)
xs = xs[::-1]
ys = ys[::-1]
#ys = ys*((xs[-1]-xs[0])/max(ys)) 

yr = max(ys)-min(ys)
xr = max(xs)-min(xs)
#ys = ys*(xr/yr) # make square (essential)

# add on dummy data beyond range to help at limits (will be removed later)
dummy_xlow = np.arange(min(xs)-rad,min(xs)-1,5)
dummy_ylow = np.zeros(np.shape(dummy_xlow)) + ys[0] *1.5
xs = np.concatenate((dummy_xlow, xs))
ys = np.concatenate((dummy_ylow, ys))

dummy_xhi = np.arange(max(xs)+1,max(xs)+rad,5)
dummy_yhi = np.zeros(np.shape(dummy_xhi)) + ys[-1] *1.5
xs = np.concatenate((xs, dummy_xhi))
ys = np.concatenate((ys, dummy_yhi))

test = False

alpha = dict()

for i in range(len(xs)):

    r = rad
    
    a = xs[i] # center x coord
    min_ix = np.where(xs == find_nearest(xs, a-r))
    min_ix = int(min_ix[0])
    max_ix = np.where(xs == find_nearest(xs, a+r))
    max_ix = int(max_ix[0])
    
    # data window between circle radius
    x = xs[min_ix+1:max_ix] # the +1 on the lower bound skips the lowest dummy dp
    y = ys[min_ix+1:max_ix] # not having +1 on the upper bound skips the highest dummy dp
                            # this ensures the bounds don't exceed the radius
                            # stopping errors later on
    
    r += 1
    
    run = True
    while(run):
        # sort window by largest intensity to speed up processing time
        unsorted = np.array((x,y))
        sort_done = np.array((unsorted[0][unsorted[1,:].argsort()], 
                              unsorted[1][unsorted[1,:].argsort()]))
        xsort = sort_done[0]
        ysort = sort_done[1]
        
        for j in range(len(y)):
            #try:
            b = ysort[j] - (r**2 - (xsort[j]-a)**2)**0.5 # center y coord
            f = b + (r**2 - (x-a)**2)**0.5 # circle y coords                
            # except FloatingPointError:
            #     test = True
            #     print(str(ysort[j] - (r**2 - (xsort[j]-a)**2)**0.5))
            #     print(str((r**2 - (x-a)**2)**0.5)
                      
            #     b = ysort[j] - (r**2 - (xsort[j]-a)**2)**0.5 # center y coord
            #     f = b + (r**2 - (x-a)**2)**0.5 # circle y coords  
                
            diff = y - f
            # if circle is below all data points
            if np.all(diff >= 0):
                alpha[i] = np.array([x, y, diff, a, b, f], dtype=object) #!!!
                run = False
                break
        
        # if no solution found, set to np.inf            
        if run == True:
            temp = np.zeros([6, len(x)]) #!!!
            temp[0,:] = temp[0,:] + x
            temp[1,:] = temp[1,:] + np.inf
            temp[2,:] = temp[2,:] + np.inf
            temp[3,:] = temp[3,:] + np.inf
            temp[4,:] = temp[4,:] + np.inf
            temp[5,:] = temp[4,:] + np.inf #!!!
            alpha[i] = temp
            run = False
            break
    

# copy data for manipulation        
beta = dict() 
beta = copy.deepcopy(alpha)  
# for each data point, choose minimum diff values for circle parameters
for k in range(len(beta)-1):
    # compares x coords
    overlap_k = np.intersect1d(beta[k][0], beta[k+1][0],
                               assume_unique=True, return_indices=True) 

    diff_k = beta[k][2][overlap_k[1]]
    diff_k2 = beta[k+1][2][overlap_k[2]]
    
    beta[k+1][2][overlap_k[2]] = np.minimum(diff_k, diff_k2)
    
    
    
t = list(beta[len(beta)-1][2]) # starting at the the RHS of the spectrum
xtemp = list(beta[len(beta)-1][0]) # for checking, will delete
ft = list(beta[len(beta)-1][5]) #!!!
t.reverse()
xtemp.reverse()
ft.reverse() #!!!
temp = list()
for p in range(len(beta)-2, -1, -1): # going through each subsequent curve for right to left
    if beta[p][0][0] != beta[p+1][0][0]:
        t.append(beta[p][2][0]) # append the left most value 
        xtemp.append(beta[p][0][0])
        ft.append(beta[p][5][0])
    else:
        temp.append(beta[p][0][0])

t = [ts + 10 for ts in t]

ft.reverse()# !!!
xtemp.reverse()
t.reverse() # reverse the spectrum

ft = np.array(ft) #!!!
xtemp = np.array(xtemp)
t = np.array(t)


# Plotting

xtemp, t = spectrum_interval(xtemp,t, lb, ub)
xs, ys = spectrum_interval(xs, ys, lb, ub)
xt, ft = spectrum_interval(xtemp, ft, lb, ub)

fig1 = plt.figure(figsize=(6.75,6.75))    
ax1 = fig1.add_subplot(111)
#label = 'r = ' + str(r), linestyle=''
ax1.plot(xs, ys, linestyle='', marker = 'o', markersize =0.3, color='r', label='Raw')
ax1.plot(xtemp, t, linestyle='', marker = 'o', markersize =0.3, color = 'k', label='RBF processed')

rbf = ys[np.intersect1d(xs, xtemp, return_indices=True)[1]] - t
ax1.plot(xtemp, rbf, linestyle='', marker = 'o', markersize =0.3, color = 'orange', label='RBF '+ str(rad))


coords = [500,800,1100,1400,1700,2000,2300,2600,2900,3200,3700]
alphas = [0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55]

alphas.reverse()


# for i in range(len(coords)):
#     coord = coords[i]
#     ac = alpha[coord][3]
#     bc = alpha[coord][4]
#     circle = plt.Circle((ac,bc), rad, color='orange', alpha=alphas[i])
#     ax1.add_patch(circle)
       

# coord = 1900
# ac = alpha[coord][3]
# bc = alpha[coord][4]
# circle = plt.Circle((ac,bc), rad, color='orange', alpha=0.3)
# ax1.add_patch(circle) 
    
# #rolling_ball check
# coord = 2500
# xc = alpha[coord][0]
# yc = alpha[coord][1]
# ac = alpha[coord][3]
# bc = alpha[coord][4]
# circle = plt.Circle((ac,bc), rad, color='orange', alpha=0.5)
# ax1.add_patch(circle)


#ax1.set_ylim(-100,1100)
ax1.legend(markerscale=8,frameon=True, loc='upper left', fancybox=False, edgecolor='k')
ax1.set_xlabel('Raman shift (cm$^{-1}$)')
ax1.set_ylabel('Intensity (arb. units)')
#ax1.set_yticks([])
#ax1.set_ylim(-100,4000)
       
# y1 = bc + (r**2 - (xc-ac)**2)**0.5
# y2 = bc - (r**2 - (xc-ac)**2)**0.5

#fig, ax = plt.subplots(figsize=(5,5))
# ax1.plot(xc, y1, color = 'k')
# ax1.plot(xc, y2, color = 'k')
# ax1.plot(xc, yc, linestyle='', marker = 'o', markersize =0.3, color = 'k')

fig1.savefig('C:/Users/'+user+'/OneDrive - University College London' + \
                    '/PhD/Python/Data/Figures/rbf_example.png', dpi=600)

fig2 = plt.figure(figsize=(6,4))
ax2 = fig2.add_subplot(111) 


xs2, t2 = spectrum_interval(xtemp,t,lb,3170)
t2 += 1000
ax2.plot(xs2, t2, linestyle='-', marker = 'o', markersize =0.8, color='k')
#ax2.set_ylim(0,40000)
#ax2.legend(markerscale=8,frameon=True, loc='upper right', fancybox=False, edgecolor='k', label='RBF processed')
ax2.set_title('Point Scan 1, RBF')
ax2.set_xlabel('Raman shift (cm$^{-1}$)')
ax2.set_ylabel('Intensity (arb. units)')
#ax2.set_yticks([])
fig2.tight_layout()
plt.axvline(x=305)
plt.axvline(x=315)
plt.axvline(x=655)
plt.axvline(x=700)

fig2.savefig('C:/Users/'+user+'/OneDrive - University College London' + \
                    '/PhD/Python/Data/Figures/rbf_processed2.svg', dpi=600, transparent = True)

plt.show()

winsound.Beep(frequency, duration)
end_time = time.perf_counter()
sys.stdout.write('\n')
print("\nScript runtime: %.2f \bs" % (end_time - start_time))

xspec = xs2
yspec = t2

