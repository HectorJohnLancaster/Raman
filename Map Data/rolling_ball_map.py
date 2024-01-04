
import numpy as np
import copy

# https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def spectrum_interval(x,y, lb, ub):
    lowerx = np.where(x>lb)
    upperx = np.where(x<ub)
    interx = np.intersect1d(lowerx,upperx)
    x = x[interx]  
    y = y[interx]
    return(x,y)

def remove_background(xs, ys, rad):

    # xs = clean_data[sample][x,y][:,2]
    # ys = clean_data[sample][x,y][:,3]  
    #xs = xsafe
    #ys = ysafe

    # xn, yn = xwindow(xs, ys, 3)
    # xs = xn
    # ys = yn
    
    lb = int(min(xs)-1) # lower and upper bounds for cut-off at end
    ub = int(max(xs)+1)

    xs = xs[::-1]
    ys = ys[::-1]

    yr = max(ys)-min(ys)
    xr = max(xs)-min(xs)
    ys = ys*(xr/yr) # make square (essential)

    # add dummy data to allow the code to catch the ends of the real data,
    # deleted later  
    dummy_xlow = np.arange(min(xs)-rad//2,min(xs)-1,5)
    dummy_ylow = np.zeros(np.shape(dummy_xlow)) + ys[0]
    xs = np.concatenate((dummy_xlow, xs))
    ys = np.concatenate((dummy_ylow, ys))
    
    dummy_xhi = np.arange(max(xs)+1,max(xs)+rad//2,5)
    dummy_yhi = np.zeros(np.shape(dummy_xhi)) + ys[-1]
    xs = np.concatenate((xs, dummy_xhi))
    ys = np.concatenate((ys, dummy_yhi))
                        
    alpha = dict()
    
    for i in range(len(xs)):
        
        r = rad
        
        a = xs[i] # center x coord
        min_ix = np.where(xs == find_nearest(xs, a-r)) # (1)
        min_ix = int(min_ix[0]) # (2)
        max_ix = np.where(xs == find_nearest(xs, a+r))
        max_ix = int(max_ix[0])
        # (1) find nearest x that corresponds to a radius r from center
        # (2) save as integer to work with slice below
        
        # data window corresponding to radius
        x = xs[min_ix+1:max_ix] # the +1 on the lower bound skips the lowest dummy dp
        y = ys[min_ix+1:max_ix] # not having +1 on the upper bound skips the highest dummy dp
                                # this ensures the bounds don't exceed the radius
                                # stopping errors later on
        
        r += 1
    
    
        run = True
        while(run):
            unsorted = np.array((x,y))
            sort_done = np.array((unsorted[0][unsorted[1,:].argsort()],
                                  unsorted[1][unsorted[1,:].argsort()]))
            xsort = sort_done[0]
            ysort = sort_done[1]
            
            for j in range(len(y)):
                b = ysort[j] - (r**2 - (xsort[j]-a)**2)**0.5 # center y coord
                f = b + (r**2 - (x-a)**2)**0.5
                
                diff = y - f
                if np.all(diff >= 0):
                    alpha[i] = np.array([x, y, diff, a, b], dtype='object')
                    run = False
                    break
                        
            if run == True:
                temp = np.zeros([5, len(x)])
                temp[0,:] = temp[0,:] + x
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
        overlap_k = np.intersect1d(beta[k][0], beta[k+1][0], 
                                   assume_unique=True, return_indices=True) # compares x coords
    
        diff_k = beta[k][2][overlap_k[1]]
        diff_k2 = beta[k+1][2][overlap_k[2]]
        
        beta[k+1][2][overlap_k[2]] = np.minimum(diff_k, diff_k2)
        
        
        
    t = list(beta[len(beta)-1][2]) # starting at the the RHS of the spectrum
    xtemp = list(beta[len(beta)-1][0]) # for checking, will delete
    t.reverse()
    xtemp.reverse()
    for p in range(len(beta)-2, -1, -1): # going through each subsequent curve from right to left
        if beta[p][0][0] != beta[p+1][0][0]:
            t.append(beta[p][2][0]) # append the left most value 
            xtemp.append(beta[p][0][0])

    

    t = [ts + 20 for ts in t]
    xtemp.reverse()
    t.reverse() # reverse the spectrum
    
    xtemp = np.array(xtemp)
    t = np.array(t)
    xtemp, t = spectrum_interval(xtemp,t, lb, ub)
    #plt.plot(xtemp, t, linestyle='', marker = 'o', markersize =0.3, color='r', label='Raw')
    #xtemp, t = spectrum_interval(xtemp,t, 1000, 1950) #!!! 
    
    return(xtemp, t) # have to add on a factor here as otherwise fit gets confused at zeros
    
    
    

