#--------------------------------File notes-----------------------------------

# Filenames must be named as "x map.txt" where x is the material ID.
# The filepath is C:\\Users\\Hector\\Desktop\\Data\\Map_Data\\ so this must
# be changed if the files to be processed are in a different location.


#---------------------------------User Inputs---------------------------------

loaded = False # load data instead of processing
name = 'CD G6 G7 785 clean data' # name of data, if any, to be loaded

# Modified z-score threshold:
threshold = 6 #6 default for Raman paper analysis

# Does the spectra have a flourescent background you wish to remove?
rad = [499]#,400]#,1000]#[1000,350] # 350 for CD7, 700 for CD6 785, 2000 for CD6 514
remove_flourescence = True
despike = True
despike_window = False
normalise = False
window = True
#combine_maps = True
b_width = 3 #3 default for Raman paper analysis

#-------------------------------Import Modules--------------------------------

import numpy as np
import time
import glob
import copy
import sys
import os
import pickle
import matplotlib.pyplot as plt

user = os.getlogin()
#------------------------------Define Equations-------------------------------

# Function to import previously stored data from pickle
def load_data(name): # name defined by user 
    path = 'C:/Users/'+user+'/OneDrive - University College London/PhD/' + \
       'Python/Data/Reports/' + name + '.pickle'  
    with open(path, 'rb') as handle:
        loaded_data = pickle.load(handle)  
        
    return loaded_data

# Equation to calculate the stepsize in physical coordinates
def stepsize(steps): # input data of coords
    trim = set(steps) # turn to a set (gets rid of mutliples)
    ltrim = list(trim) # turn into a list (gives positional meaning)
    ltrim.sort() # sorts the list
    if len(ltrim) < 2:
        ltrim.append(ltrim[0])
    return ltrim[1] - ltrim[0] # return a step size

# Extract information about the map setup
def mapinfo(sample):
    xmin = min(raw_data[sample][:,1]) # see data and image for column index
    xmax = max(raw_data[sample][:,1])
    ymin = min(raw_data[sample][:,0])
    ymax = max(raw_data[sample][:,0])
    xstep = stepsize(raw_data[sample][:,1]) # finds the step size
    ystep = stepsize(raw_data[sample][:,0])
    return xmin, xmax, ymin, ymax, xstep, ystep
    
import winsound
frequency = 500  # Set Frequency To 2500 Hertz
duration = 1000  # Set Duration To 1000 ms == 1 second


#--------------------------------Import Data----------------------------------

# Print code status info
sys.stdout.write('\nImporting data: ...')

# Start process timer
start_time = time.perf_counter() 

# Create a list of all files in filepath that have the form *.txt
filenames = sorted(glob.glob(r'C:/Users/'+user+'/OneDrive - University College' + \
                   ' London/PhD/Python/Data/Map_Data/*.txt'))
    
idx = filenames[0].find('\\') + 1

# Extract the sample names from the filepath name
sample_names = list()
for files in range(len(filenames)):
    sample_names.append(filenames[files][idx:-4]) # Finds position of name


# This numpy method loads the text in as an array, which is matched to the 
# corresponding material key and stored in the dictionary 'raw_data'   
raw_data = dict()
for sample in range(len(filenames)):
        raw_data[sample_names[sample]] = np.loadtxt(filenames[sample])

# Print status info
sys.stdout.write('\rImporting data: 100%\n')


#----------------------------------Load Data----------------------------------

if loaded == True:
    sys.stdout.write('\nLoading data: ...')
    #clean_data = load_data(name)
    
    #---increase bin size---
    if window == True:
        from xwindow import xwindow # increases bin size by averaging data over i+-b_width
        for sample in sample_names:
            xmin, xmax, ymin, ymax, xstep, ystep = mapinfo(sample)  
            for x in np.arange(xmin, xmax+xstep, xstep):
                x = np.round(x, decimals = 2)
                for y in np.arange(ymin, ymax+ystep, ystep):
                    y = np.round(y, decimals = 2)
                    xs = clean_data[sample][(x,y)][2,:]
                    ys = clean_data[sample][(x,y)][3,:]                    
                    if np.isnan(ys).all() != True:
                        xbin, ybin = xwindow(xs,ys, b_width)
                        length = len(xbin)
                        coordy = np.zeros(length) + y
                        coordx = np.zeros(length) + x  
                        coord = (x,y)
                        new_data = np.array([coordx, coordy, xbin, ybin])
                        clean_data[sample][coord] = new_data
                        
    sys.stdout.write('\rLoading data: 100%\n')
    sample_names = list(clean_data.keys())
else:    
#---------------------------------Slice Data----------------------------------

    s_count = 0
    sample_num = str(len(sample_names))
    sliced_data = dict()
    for sample in sample_names:
        xmin, xmax, ymin, ymax, xstep, ystep = mapinfo(sample)
        spectra = dict()
        # --- percentage counter bits ---
        pc_count = 0 # initialise
        pc_tot = len(np.arange(xmin, xmax+xstep, xstep)) # denominator in % calc
        s_count += 1 # counts the number of samples
        sys.stdout.write('\n')
        # -------------------------------
        for x in np.arange(xmin, xmax+xstep, xstep): # for each xcoord
            x = np.round(x, decimals = 2)
            # --- percentage counter bits ---
            pc_count += 1 # progress counter for % calc
            pc = (pc_count*100)//(pc_tot) # % calc
            sys.stdout.write('\rProcessing ' + str(s_count) + ' of ' + \
                             sample_num + ': ' + str(pc) + '%')
            sys.stdout.flush()
            # -------------------------------
            for y in np.arange(ymin, ymax+ystep, ystep): # for each ycoord
                y = np.round(y, decimals=2)
                posx = np.where(raw_data[sample][:,1] == x) # coords when true
                posy = np.where(raw_data[sample][:,0] == y) # coords when true
                inter = np.intersect1d(posx,posy) # intersection of posx & posy
                #-------------------------------------------------------------
                # only for when the txt file is a mosiac of overlapping maps
                idxs = np.where(np.diff(inter)!=1)[0] 
                if len(idxs) != 0:
                    idx = np.where(np.diff(inter)!=1)[0][0] + 1
                    inter = inter[0:idx]
                #-------------------------------------------------------------   
                spectra[x,y] = raw_data[sample][inter] # assign spectrum to dict
        sliced_data[sample] = spectra # for the sample, embed spectra dict in dict
    
    
    #--------------------------------Data Cleaning----------------------------
    
    # Method adapted from: 
    # Whitaker, D and Hayes, K. “A simple algorithm for despiking Raman spectra.”
    # Chemometrics and Intelligent Laboratory Systems 179 (2018): 82–84.
    
    
    #---define zscore---
    
    def modified_z_score(xs, ys, lowsafe, highsafe): # find the modified z-score of a dataset
        ydif_hl = np.diff(ys) # (1)
        ydif_hl = np.append(ydif_hl, 0)
        ydif_lh = np.diff(ys) # (2)
        ydif_lh = np.append(0, ydif_lh)
        ydif = np.maximum(np.abs(ydif_hl), np.abs(ydif_lh)) # (3)
        median = np.median(ydif) # median of these differences
        mad = np.median([np.abs(ydif-median)]) # median absolute deviation (MAD)
        if mad == 0:
            z = np.empty(np.shape(ydif))
            z[:] = np.nan
            plt.plot(xs, ys)
            plt.title('clean data set to nan')
        else:
            z = 0.6745 * (ydif - median) / mad # the z-score, 0.6745 adjusting factor
            
        safe = np.logical_and(xs>lowsafe, xs<highsafe)==False
        z = z*safe + 1  # if zscore is within safe region, will be *False = 0, 
                        # +1 means all safe values will have a zscore of 1
        
        return z
    
    # (0) General note: the way the code calculates the differences in f has been
    #     changed from the referenced link. I believe this way solves an issue
    #     whereby peaks with similar intensities spanning multiple points were not
    #     always flagged. Before the code only calculated differences as is 
    #     described in (1) below. The issue now seems to be resolved.
    # (1) The high to low difference 'hl' in f. A zero added on at the end of this
    #     list as the last number has no higher number to compare difference to.
    # (2) The low to high differences 'lh' in f. The first value is then set to 
    #     zero as there is no number below the first entry (and the way python
    #     indexes will subtract the last number in f, which we don't want)
    # (3) Compare the hl and lh differences, choose the largest number.
    
    
    
    fail_info = set()
    #---define fixer---
    def fixer(f, z_bool, ma=10): # ma sets the 'data window' size
        n = len(z) # for indexing the arrays
        f_out = f.copy()
        z_bool[0] = z_bool[n-1] = True # set the first and last value to True
        spikes = np.array(np.where(z_bool == True)) # gives spike locations
        spikes = np.reshape(spikes,(len(spikes[0]),)) # (1)
        for j in spikes: # uses the locations of spikes as j values
            w = np.arange(max(0, j - ma), min(n, j + 1 + ma)) # (2)
            w = w[z_bool[w] == 0] # returns the non-peak locations from range
            if len(w) != 0:
                f_out[j] = np.mean(f[w]) # (3)
            else:
                global fail_info
                fail_info.add(sample +': ' + str((x,y)))
                
        return [f_out, len(spikes)] 
    
    # (1) must reshape from (1,len) to (len,) for the cycling through spike
    # location 'j' in the loop to work.
    # (2) sets a range 'w' around the spike location 'j', if close to array
    # boundaries, limits range to boundary, +1 becuase way python references.
    # (3) replace peak with mean of surrounding non-peak data.
    
    
    # #---clean data---
    
    if despike == True:
        # num_spikes = dict() this would be to count number of spikes removed
        #clean_data = copy.deepcopy(sliced_data) # otherwise overwrites
        clean_data = dict()
        spike_count = 0
        if despike_window == True:
            lowsafe = float(input("Enter lower wavenumber limit for this region: "))
            highsafe = float(input("Enter upeer wavenumber limit for this region: "))
        else:
            lowsafe = 0
            highsafe = 1
        for sample in sample_names:
            xmin, xmax, ymin, ymax, xstep, ystep = mapinfo(sample)
            spectra = dict()
            for x in np.arange(xmin, xmax+xstep, xstep):
                x = np.round(x, decimals = 2)
                for y in np.arange(ymin, ymax+ystep, ystep):
                    y = np.round(y, decimals = 2)
                    
                    ys = sliced_data[sample][(x,y)][:,3] #(1)                    
                    xs = sliced_data[sample][(x,y)][:,2]
                    z = modified_z_score(xs, ys, lowsafe, highsafe) # (2) 
                    z_bool = np.abs(z) > threshold # (3)
                    
                    if np.isnan(z).all() != True:
                        ys_fixed, spikes = fixer(ys,z_bool) # send to fixer function
                        #ys_fixed += 10000 # helps with fitting
                        spike_count += spikes
                        length = len(ys_fixed)
                        coordy = np.zeros(length) + y
                        coordx = np.zeros(length) + x      
                    else:
                        length = len(ys)
                        coordy = np.zeros(length) + y
                        coordx = np.zeros(length) + x
                        ys[:] = np.nan
                        ys_fixed = ys
                        xs[:] = np.nan
                        
                    spectra[x,y] = np.array([coordx, coordy, xs, ys_fixed])
                    
            clean_data[sample] = spectra
            
    else:
        # redo entry just to allow df to be in consistant format
        clean_data = dict()
        for sample in sample_names:
            xmin, xmax, ymin, ymax, xstep, ystep = mapinfo(sample)
            spectra = dict()
            for x in np.arange(xmin, xmax+xstep, xstep):
                x = np.round(x, decimals = 2)
                for y in np.arange(ymin, ymax+ystep, ystep):
                    y = np.round(y, decimals = 2)
                    
                    ys = sliced_data[sample][(x,y)][:,3]# + 10000 # helps with fitting           
                    xs = sliced_data[sample][(x,y)][:,2]
                    coordx = sliced_data[sample][(x,y)][:,1]
                    coordy = sliced_data[sample][(x,y)][:,0]
                    spectra[x,y] = np.array([coordx, coordy, xs, ys])
        
            clean_data[sample] = spectra
    
    
    # # (1) create (n,) array of intensity values
    # # (2) send f to zscore function
    # # (3) create (n-1,) array of booleans which give true if spike detected
    # #     in the data, correspinding to a zscore greater than the threshold.
    # # (4) update column w/ despiked data
    
    
    #---increase bin size---
    if window == True:
        from xwindow import xwindow # increases bin size by averaging data over i+-b_width
        for sample in sample_names:
            xmin, xmax, ymin, ymax, xstep, ystep = mapinfo(sample)  
            for x in np.arange(xmin, xmax+xstep, xstep):
                x = np.round(x, decimals = 2)
                for y in np.arange(ymin, ymax+ystep, ystep):
                    y = np.round(y, decimals = 2)
                    xs = clean_data[sample][(x,y)][2,:]
                    ys = clean_data[sample][(x,y)][3,:]                    
                    if np.isnan(ys).all() != True:
                        xbin, ybin = xwindow(xs,ys, b_width)
                        length = len(xbin)
                        coordy = np.zeros(length) + y
                        coordx = np.zeros(length) + x  
                        new_data = np.array([coordx, coordy, xbin, ybin])
                        coord = (x,y)
                        clean_data[sample][coord] = new_data
                        
                    
        
    #------------------------------Normalisation-------------------------------
    
    if normalise == True:
        norm_data = copy.deepcopy(clean_data)
        for sample in sample_names:
            xmin, xmax, ymin, ymax, xstep, ystep = mapinfo(sample)
            spectra = dict()
            for x in np.arange(xmin, xmax+xstep, xstep):
                x = np.round(x, decimals = 2)
                for y in np.arange(ymin, ymax+ystep, ystep):
                    y = np.round(y, decimals = 2)
                    f = norm_data[sample][(x,y)][3,:] # intensity data column
                    xs = norm_data[sample][(x,y)][2,:] # wavenumber data column
                    
                    #===Remove a linear background before normalization===
                    # m = (f[-1]-f[0])/(xs[-1]-xs[0])
                    # c = f[-1]-m*xs[-1]
                    # linback = m*xs + c
                    # f = f - linback
                    #=====================================================
                    
                    lowerx = np.where(xs>1250) # define region to normalise
                    upperx = np.where(xs<1350) # with lower and upper bounds
                    interx = np.intersect1d(lowerx,upperx) # find xlocs in region      
                    max_p = max(f[interx]) # maximum intensity in region
                    f_norm = f/max_p # normalise to max_p w/ value = 1
                    norm_data[sample][(x,y)][3,:] = f_norm # replce intensity data
                                
    
    #----------------------------Flourescence Removal-------------------------
    
    
    if remove_flourescence == True:
        counter = 0
        sys.stdout.write('\n')
        s_count = 0
        import rolling_ball_map
        
        for sample in sample_names:
    
            xmin, xmax, ymin, ymax, xstep, ystep = mapinfo(sample)
            
            # --- percentage counter bits ---
            pc_count = 0 # initialise
            lenx = len(np.arange(xmin, xmax+xstep, xstep))
            leny = len(np.arange(ymin, ymax+ystep, ystep))
            pc_tot = lenx*leny # denominator in % calc
            s_count += 1 # counts the number of samples
            sys.stdout.write('\n')
            # -------------------------------
            
            r = rad[counter]
            print('radius = ' + str(r))
            
            for x in np.arange(xmin, xmax+xstep, xstep):
                x = np.round(x, decimals = 2)         
                for y in np.arange(ymin, ymax+ystep, ystep):
                    # --- percentage counter bits ---
                    pc_count += 1 # progress counter for % calc
                    pc = (pc_count*100)//(pc_tot) # % calc
                    sys.stdout.write('\rRemoving flourescent background ' + \
                                     str(s_count) + ' of ' + sample_num + \
                                         ': ' + str(pc) + '%')
                    sys.stdout.flush()
                    # ------------------------------- 
                    
                    y = np.round(y, decimals = 2)
                    
                    # xs = clean_data[sample][x,y][:,2]
                    # ys = clean_data[sample][x,y][:,3]
                    xs = clean_data[sample][x,y][2,:]
                    ys = clean_data[sample][x,y][3,:] 
                    
                    if np.isnan(ys).all() != True:
                        
                        x_mod, y_mod = rolling_ball_map.remove_background(xs,
                                                                          ys,
                                                                          r)                    

                        #y_mod = y_mod + 10000 
                        # as otherwise the baseline at zero has trouble when fitting
                        length = len(x_mod)
                        coordy = np.zeros(length) + y
                        coordx = np.zeros(length) + x
                        coord = (x,y)
                        new_data = np.array([coordx, coordy, x_mod, y_mod])
                
                        clean_data[sample][coord] = new_data
                        
            counter += 1


# if combine_maps == True:
#     for sample in sample_names:
#         if sample == sample_names[0]:
#             dc = clean_data[sample].copy() # combined dict
#             print(len(dc))
#         else:
#             dx = clean_data[sample] # dict x 
#             dc.update(dx) # add dict x to combined dict
#             print(len(dx))
#             print(len(dc))

        

                
#-----------------------------------------------------------------------------

# End process timer
end_time = time.perf_counter()
sys.stdout.write('\n')
print("\nScript runtime: %.2f \bs" % (end_time - start_time))
winsound.Beep(frequency, duration)
# last runtime = 120s

#---------------------------------Script End----------------------------------
