#--------------------------------File notes-----------------------------------

# Filenames must be named as "x map.txt" where x is the sample_names ID.
# The filepath is C:\\Users\\Hector\\Desktop\\Data\\Map_Data\\ so this must
# be changed if the file is in a different location.

# This code produces a library of spectra with the keys being the materials'
# names, and within those keys a second key being the xy position. e.g. to 
# access KC10 spectra x=0, y=0 then type: sliced_data['sample_names ID'][0,0]

#---------------------------------User Inputs---------------------------------

# Modified z-score threshold:
threshold = 1

# Enter sample_names name for save files
#CB 1400
#CD2 2000 x2
#CD6 1000
#CD7 350
#CD8 350 x2
#Graphite 1000
#YP50 2000
mat_nam = 'raw'
#rad = [2000, 350, 2000, 2000]#, 2000] # 700 not 1000 for CD6
rad = [400]#,1200,1200]#, 100, 2000, 2000, 700, 350, 350, 350, 1000, 2000]
despike = False
window = True
remove_flourescence = False
b_width = 4



#-------------------------------Import Modules--------------------------------
from package.Raman import * # imports everything in the Raman package 
                            # this includes various modules and functions 
import numpy as np
import glob
import copy
import matplotlib.pyplot as plt
import os

user = os.getlogin()

def spectrum_interval(x,y, lb, ub):
    lowerx = np.where(x>lb)
    upperx = np.where(x<ub)
    interx = np.intersect1d(lowerx,upperx)
    x = x[interx]  
    y = y[interx]
    return(x,y)

def normalise(xs, ys):
    maxi = max(ys)
    m = (np.mean(ys[-10:] - np.mean(ys[:10])))/(xs[-1]-xs[0])
    c = ys[-1]-m*xs[-1]
    linback = m*xs + c
    ys = (ys-linback)/maxi
    return(xs,ys)

import winsound
frequency = 500  # Set Frequency To 2500 Hertz
duration = 1000  # Set Duration To 1000 ms == 1 second

#--------------------------------Import Data----------------------------------

# Create a list of all files in filepath that have the form *.txt
filenames = sorted(glob.glob('C:/Users/'+user+'/OneDrive - University' + \
                             ' College London/PhD/Python/Data/Point_Data/*.txt'),
                   reverse = False)

# Extract the sample_names names from the filepath name
sample_names = list()
for sample in range(len(filenames)):
    sample_names.append(filenames[sample][80:-4]) # Corresponds to position of name

# This numpy method loads the text in as an array   
raw_data = dict()
for sample in range(len(filenames)):
        raw_data[sample_names[sample]] = np.loadtxt(filenames[sample])


#--------------------------------Data Cleaning--------------------------------

# Method adapted from: 
# Whitaker, D and Hayes, K. “A simple algorithm for despiking Raman spectra.”
# Chemometrics and Intelligent Laboratory Systems 179 (2018): 82–84.


#---define zscore---

def modified_z_score(f): # Find the modified z-score of a dataset f
    median = np.median(f)
    mad = np.median([np.abs(f-median)]) # Median absolute deviation (MAD)
    z = 0.6745 * (f - median) / mad # The z-score, 0.6745 adjusting factor
    return z


#---define fixer---

def fixer(f, z_bool, ma=10):
    n = len(z) # for indexing the arrays, (575,)
    f_out = f.copy()
    z_bool[0] = z_bool[n-1] = True # set the first and last value to True
    spikes = np.array(np.where(z_bool == True)) # gives spike locations
    spikes = np.reshape(spikes,(len(spikes[0]),)) # (1)
    for j in spikes: # uses the locations of spikes as j values
        w = np.arange(max(0, j - ma), min(n, j + 1 + ma)) # (2)
        w = w[z_bool[w] == 0] # returns the non-peak locations from the range
        f_out[j] = np.mean(f[w]) # (3)
        #f_out[j] = np.nan
    return f_out # (576,)

# (1) must reshape from (1,len) to (len,) for the cycling through j in the for
# loop to work.
# (2) sets a range around j, if close to array boundaries, limits range to
# boundary, +1 becuase the way python references does not include top value.
# (3) replace peak with mean of surrounding non-peak data.


#---clean data---

# num_spikes = dict() this would be to count number of spikes/extremes removed
if despike == True:
    clean_data = copy.deepcopy(raw_data) # need deepcopy otherwise overwrites
    for sample in sample_names:
        spectra = dict()
        f = raw_data[sample][:,1] # (1)
        z = modified_z_score(np.diff(f)) # (2)   
        z_bool = np.abs(z) > threshold # (3)
        f_fixed = fixer(f,z_bool) # send to fixer function
        clean_data[sample][:,1] = f_fixed # (4)
else:
    clean_data = copy.deepcopy(raw_data)  

# (1) create a (n,) array of intensity values f
# (2) send a (n-1,) array of differences in f to zscore function
# (3) create a (n-1,) array of booleans which give true if a spike is detected
#     in the data, correspinding to a zscore greater than the threshold.
# (4) update column w/ despiked data

#---increase bin size---
if window == True:
    from xwindow import xwindow # increases bin size by averaging data over i+-b_width
    for sample in sample_names:    
        x = clean_data[sample][:,0]
        y = clean_data[sample][:,1]                    
        if np.isnan(y).all() != True:
            xbin, ybin = xwindow(x,y, b_width)
            length = len(xbin)
            new_data = np.zeros((length,2))
            
            new_data[:,0] = xbin
            new_data[:,1] = ybin
            #new_data = np.array([xbin, ybin])
            clean_data[sample] = new_data     

        
#--------------------------------Normalisation---------------------------------

norm_data = dict()
for sample in sample_names:
    x = clean_data[sample][:,0]           
    y = clean_data[sample][:,1]
    x, y = spectrum_interval(x,y,150,3200) # set normalise range
    x, f_norm = normalise(x, y)
    #f_norm = y/max(f) # normalise to highest point in range
    norm_data[sample] = np.zeros((len(x), 2))
    norm_data[sample][:,1] = f_norm
    norm_data[sample][:,0] = x
            
    
#----------------------------Flourescence Removal-----------------------------

counter=0
if remove_flourescence == True:
    s_count = 0
    import rolling_ball_single
    for sample in sample_names:
        x = clean_data[sample][:,0]
        y = clean_data[sample][:,1]
        r = rad[counter]
        x_mod, y_mod = rolling_ball_single.remove_background(x, y, r)
        
        length = len(x_mod)
        new_data = np.array([x_mod, y_mod])
        
        clean_data[sample] = new_data   

        counter += 1         

    # x_mod, y_mod = spectrum_interval(x_mod,y_mod,1000,1950)
    # plt.plot(x_mod, y_mod)

fig, ax = plt.subplots(figsize=(4,4))
counter = 0
#E69F00
#56B4E9
#009E73
#F0E442
#CC79A7
c2 = '#56B4E9'
c3 = '#009E73'
c4 = '#F0E442'
c5 = '#CC79A7'
c1 = '#E69F00'
colors = [c2,c1,c3,c4,c5]
gap = 0
for sample in sample_names:            
    # x = norm_data[sample][:,0]           
    # y = norm_data[sample][:,1]
    x = norm_data[sample][:,0]           
    y = norm_data[sample][:,1]
    x, y = spectrum_interval(x,y+gap,150,3200)
    #counter += 800
    plt.plot(x,y, linestyle='-', marker='', markersize=1.2, mew=0, 
             label=sample, color=colors[counter])
    counter+=1
    gap+=1

#plt.legend()
plt.yticks([])
plt.xlim(1,3350)
#plt.axvline(x=1450, linestyle = ':', marker = '', color = 'k', linewidth = 0.9, zorder=20)
ax.xaxis.set_major_locator(ticker.MultipleLocator(500))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(100))
plt.xlabel('Raman shift (cm$^{-1}$)')
plt.ylabel('Intensity (arb. units)')
xsafe = x #different ids for rb checker
ysafe = y
plt.tight_layout(pad=1.08, h_pad=None, w_pad=None, rect=None)

plt.savefig('C:/Users/'+user+'/OneDrive - University College London' + \
            '/PhD/Python/Data/Figures/'+' '.join(sample_names)+' fullscan stacked.tif', transparent=True, dpi=300)
plt.savefig('C:/Users/'+user+'/OneDrive - University College London' + \
            '/PhD/Python/Data/Figures/'+' '.join(sample_names)+' fullscan stacked.svg', transparent=True, dpi=600)
    
#winsound.Beep(frequency, duration)

