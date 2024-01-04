
#--------------------------------file notes-----------------------------------

# Creates a combined plot containing all map g peak postition wrt their xy 
# coordinates for each material. 
    

#------------------------------import modules---------------------------------

import matplotlib as mpl
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec
import copy
import time
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import FixedLocator, FixedFormatter

def spectrum_interval(x,y, lb, ub):
    lowerx = np.where(x>lb)
    upperx = np.where(x<ub)
    interx = np.intersect1d(lowerx,upperx)
    x = x[interx]  
    y = y[interx]
    return(x,y)

def mapinfo(sample):
    xmin = min(raw_data[sample][:,1]) # see data and image for column index
    xmax = max(raw_data[sample][:,1])
    ymin = min(raw_data[sample][:,0])
    ymax = max(raw_data[sample][:,0])
    xstep = stepsize(raw_data[sample][:,1]) # finds the step size
    ystep = stepsize(raw_data[sample][:,0])
    return xmin, xmax, ymin, ymax, xstep, ystep

def stepsize(steps): # input data of coords
    trim = set(steps) # turn to a set (gets rid of mutliples)
    ltrim = list(trim) # turn into a list (gives positional meaning)
    ltrim.sort() # sorts the list
    if len(ltrim) < 2:
        ltrim.append(ltrim[0])
    return ltrim[1] - ltrim[0] # return a step size

# from contextlib import contextmanager
# @contextmanager
# def temp_rcParams(params):
#     backup = {key:plt.rcParams[key] for key in params}  # Backup old rcParams
#     plt.rcParams.update(params) # Change rcParams as desired by user
#     try:
#         yield None
#     finally:
#         plt.rcParams.update(backup) # Restore previous rcParams
        

#--------------------------------user inputs----------------------------------

thresh_Chi = np.inf#10#0.037#0.13#0.037#0.03#0.037#0.04#0.038#35 # Can change here to tune
thresh_sig = np.inf#50#2 #0.05 #5 is default minimum cutoff for thesis (don't go below)
med_num = 5 # If

#title = 'PosG'

#sample_names = ['Char kc10 map']#, 'CB KC24 map', 'CB KC24 quenched map']

savefigures = True
median_filter = False 
hard_lim = False
log = False

if hard_lim == True:
    limmax = 500#np.inf#1605#1525#1
    limmin = 0.0001#-np.inf#1545#1440 #1
else:
    limmax = np.inf
    limmin = -np.inf

nan_var = 'cen2'

# choose the fit variable to map
fit_var = ['Area 2ovr3']#,'cen2','cen3','cen4','Area 4ovr2']#, 'cen2', 'Gamma2']

# interpolation: choose either "lanczos" or "none"
inter = "none"

# format scalebar
scalebar_microns = 20 
scalebar_height = 0.03
textsize = 10

# choose if figure is half or full page width
full_width = 6.5
half_width = 3.25
width = 6.5#*1.095967684107369

# scale the overall image in its container wrt full_width
size_scale = 1#0.9#0.5

# v=verticle, h=horizontal
orientation = 'h' 

#----------------------------start process timer------------------------------

start_time = time.perf_counter()


#--------------------------format layout variables----------------------------

pc_scale = 0.85
if orientation == 'v':
    hspace = 0.1
    wspace = 0.0
    top = 0.95
    right = 0.95
    bottom = 0.15
    left = 0.15
    
    nrows = len(sample_names)+1
    ncols = 1
    height = width*size_scale
    colcount = 1
    rowcount = len(sample_names)
    
if orientation == 'h':
    hspace = 0.0
    wspace = 0.1
    top = pc_scale
    right = pc_scale
    bottom = 0.0
    left = 0.0

    nrows = 1
    ncols = len(sample_names)+1
    width = width*size_scale
    colcount = len(sample_names)
    rowcount = 1


if len(sample_names) == 1:   
    if orientation == 'h':
        hr = [1]
        wr = [1,0.08]
        widthspace = wspace*(sum(wr)/(ncols+1))
        height = 1/(sum(wr)+ncols*widthspace) * width
        
    elif orientation == 'v':
        hr = [1, 0.08]
        wr = [1]
        heightspace = hspace*(sum(hr)/nrows)
        width = 1/(sum(hr)+(nrows-1)*heightspace) * height   
    
    
elif len(sample_names) == 2:    
    if orientation == 'h':
        hr = [1]
        wr = [1,1,0.08]
        widthspace = wspace*(sum(wr)/(ncols+1))
        height = 1/(sum(wr)+ncols*widthspace) * width
        
    elif orientation == 'v':
        hr = [1,1,0.08]
        wr = [1]
        heightspace = hspace*(sum(hr)/nrows)
        width = 1/(sum(hr)+(nrows-1)*heightspace) * height    
        

elif len(sample_names) == 3:
    if orientation == 'h':
        hr = [1]
        wr = [1,1,1,0.08]
        widthspace = wspace*(sum(wr)/(ncols+1))
        height = 1/(sum(wr)+ncols*widthspace) * width
        
    elif orientation == 'v':
        hr = [1,1,1,0.08]
        wr = [1]
        heightspace = hspace*(sum(hr)/nrows)
        width = 1/(sum(hr)+(nrows-1)*heightspace) * height   
        
elif len(sample_names) == 4:    
    if orientation == 'h':
        hr = [1]
        wr = [1,1,1,1,0.08]
        widthspace = wspace*(sum(wr)/(ncols+1))
        height = 1/(sum(wr)+ncols*widthspace) * width
        
    elif orientation == 'v':
        hr = [1,1,1,1,0.08]
        wr = [1]
        heightspace = hspace*(sum(hr)/nrows)
        width = 1/(sum(hr)+(nrows-1)*heightspace) * height   

elif len(sample_names) == 5:    
    if orientation == 'h':
        hr = [1]
        wr = [1,1,1,1,1,0.08]
        widthspace = wspace*(sum(wr)/(ncols+1))
        height = 1/(sum(wr)+ncols*widthspace) * width
        
    elif orientation == 'v':
        hr = [1,1,1,1,1,0.08]
        wr = [1]
        heightspace = hspace*(sum(hr)/nrows)
        width = 1/(sum(hr)+(nrows-1)*heightspace) * height  
        
elif len(sample_names) == 7:    
    if orientation == 'h':
        hr = [1]
        wr = [1,1,1,1,1,1,1,0.08]
        widthspace = wspace*(sum(wr)/(ncols+1))
        height = 1/(sum(wr)+ncols*widthspace) * width
          

elif len(sample_names) == 13:    
    if orientation == 'h':
        hr = [1]
        wr = [1,1,1,1,1,1,1,1,1,1,1,1,1,0.08]
        widthspace = wspace*(sum(wr)/(ncols+1))
        height = 1/(sum(wr)+ncols*widthspace) * width
        
elif len(sample_names) == 15:    
    if orientation == 'h':
        hr = [1]
        wr = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.08]
        widthspace = wspace*(sum(wr)/(ncols+1))
        height = 1/(sum(wr)+ncols*widthspace) * width
        

#----------------------------scalebar properties------------------------------

bar_length = list()
bar_height = list()
for sample in sample_names:        
    xmin, xmax, ymin, ymax, xstep, ystep = mapinfo(sample) 
    square = 1/(1 + (xmax-xmin)/xstep) # size of a square
    # scalebar properties
    bar_length.append((100*square/xstep)*(scalebar_microns/100))
    bar_height.append(scalebar_height)#0.03
    

#-----------------------------process fit data--------------------------------

plt.rcParams.update({'figure.max_open_warning': 0}) # supresses runtime info

mat_grid = dict() # initialise variable
chuck = dict()
if fit_var == [nan_var]:
    nan_dict = dict()

chucker = 0


for variable in fit_var:
    chucker += 1
    for sample in sample_names:
        
        if fit_var == [nan_var]:
            nan_dict[sample] = list()
            
        if sample == sample_names[0]:#or sample == sample_names[1]:
            variable = 'Area 1ovr2'
        elif sample == sample_names[1]:
            variable = 'Area 1ovr32'
        else:
            variable = 'Area 2ovr3'
            
        # if sample == sample_names[1] or sample == sample_names[2]:
        #     variable = 'cen6'
        # else:
        #     variable = fit_var[0]
            
        xmin, xmax, ymin, ymax, xstep, ystep = mapinfo(sample) 
        # xmin = -22
        # xmax = 20   
        if sample == 'top1 doped small map1':
            xmin = -10
            xmax = 10
            
    
        #---creat dataframe of coordinates---
        temp = fit_data[sample]   
        coords = list(temp.columns) # list of all coords (x,y)
        xcoords = list()
        ycoords = list()
        for j in range(len(coords)):
            xcoords.append(coords[j][0]) # take all x coords
            ycoords.append(coords[j][1])# take all y coords
        xset = list(set(xcoords)) # make set to get rid of multiples
        yset = list(set(ycoords)) # then list so that it can be sorted below
        xset.sort() # sort in ascending order
        yset.sort()
        grid = pd.DataFrame(index = yset, columns = xset, dtype='float64')#(2)
        
        # (1) make a temporary dataframe of all fit data without variable 
        #     names, this is to allow for data manipulation (removes string)
        # (2) make a 2D grid of x and y coordinates
        #------------------------------------
        
        #---get relevant fit data---
        df = fit_data[sample] # for ease of writing
        data = df.loc[variable] # (1)
        
        # (1) make dataframe (data) with only the fit variable of interest
        # (2) drop the parameters column to allow for data manipulation
        #---------------------------
        
        #---tweak fit data---
        if variable == 'q':
            data = abs(1/data) # (1)
            axis_label = 'abs. Fano factor $| 1/q |$'
            decimals = '%.3f'
        elif variable[0:5] == 'Gamma':
            data = data 
            axis_label = 'FWHM (cm$^{-1}$)'
            decimals = '%.0f'
        elif variable[0:2] == 'I0':
            data = data
            axis_label = 'Intensity (counts/second)'  
        elif variable[:4] == 'Area':
            data = data
            p1 = str(variable[5])
            p2 = str(variable[9])
            axis_label = 'A$_'+p1+'$/A$_'+p2+'$'
            axis_label = 'A$_D$/A$_G$'
            #axis_label = 'A$_{LFP}$/A$_G$'
            decimals = '%.0f'
        elif variable[:9] == 'Intensity':
            data = data
            p1 = str(variable[-5])
            p2 = str(variable[-1])
            axis_label = 'I$_'+p1+'$/I$_'+p2+'$'
            axis_label = 'I$_D$/I$_G$'
            #axis_label = 'A$_{Mn}$$^{3+}$/A$_{Mn}$$^{4+}$'
        else:
            data = data
            axis_label = 'Raman shift (cm$^{-1}$)'
            #axis_label = 'PosD (cm$^{-1}$)'
            decimals = '%.0f'
            
        #axis_label = 'A$_{D\'\'}$/A$_G$'
        #axis_label = 'I$_{G}$/I$_{diamond}$'
        #axis_label = 'I$_{G}$/I$_{LFP}$'
        # (1) this is the fano factor, 1/q, take the absolute value
        #--------------------     
    
        #---reduced chi2---
        chi = df.loc['red_Chi2']
        chi_na = chi[data < 1E500]
        chi_na = chi_na.dropna(axis=0)
        #------------------
        
        #---data stats---
        dstat = data[data < 1E500]
        dstat = dstat.dropna(axis=0) # drop any NaNs from failed fits
        dstat = np.array(dstat)[np.where(chi_na<thresh_Chi)] # ignore'bad'data
        x0, sigma = stats.norm.fit(dstat) 
        med = np.median(dstat)
        thresh_med = med_num/med
        #----------------
        
        #---populate grid---
        chuck[sample] = list()

        for x in np.arange(xmin, xmax+xstep, xstep):
            x = np.round(x, decimals = 2)
            for y in np.arange(ymin, ymax+ystep, ystep):
                y = np.round(y, decimals = 2)
                # if y == 4:
                #     x += 2
                dp = float(data[x,y])
                dpchi = float(chi[x,y])
                dpdata = float(data[x,y])
                if dpchi > thresh_Chi or dpchi < 0.0:#0.85: # (1)
                    print('nan1')
                    dp = np.nan
                    chuck[sample].append((x,y)) # (2)
                    if variable == nan_var:
                        nan_dict[sample].append((x,y))
                if dpdata != np.nan:
                    # print('t')
                    if dpdata > x0 + thresh_sig*sigma or \
                        dpdata < x0 - thresh_sig*sigma: # (3)
                        print('nan2.1')
                        dp = np.nan
                        chuck[sample].append((x,y))
                        if variable == nan_var:
                            print('nan2')
                            nan_dict[sample].append((x,y))                    
                    if dpdata > limmax or dpdata < limmin: # (3)
                        print('nan3.1')
                        dp = np.nan
                        chuck[sample].append((x,y))
                        if variable == nan_var:
                            print('nan3')
                            nan_dict[sample].append((x,y))                    
                    if median_filter == True:
                        if dpdata < med - thresh_med or dpdata > med + thresh_med:                        
                            dp = np.nan
                            chuck[sample].append((x,y))                        
                            if variable == nan_var:
                                print('nan4')
                                nan_dict[sample].append((x,y))                        

                    
                grid[x][y] = dp
        
        #if variable == 'cen1' or variable == 'cen2' or variable == 'cen3' or variable == 'Area 4ovr2':        
        for c in range(len(nan_dict[sample])):
            print('t5' + str(x) + str(y))
            x = nan_dict[sample][c][0]
            y = nan_dict[sample][c][1]
            grid[x][y] = np.nan

        
        # (1) if the datapoint (dp) of interest lies outwith the reduced Chi2
        #     threshold, then set to NaN
        # (2) record any spectra that have been chucked, for visual inspection
        # (3) if the dp lies outwith the threshold number of standard 
        #     deviations from the data set left by (1), then set to NaN
    
        
        mat_grid[sample] = grid # assign grid to material
    

            
    #------------------------------find range---------------------------------
    
    maxlist = list()
    minlist = list()
    for sample in sample_names:
        minlist.append(np.nanmin(np.array(mat_grid[sample]))) # (1)
        maxlist.append(np.nanmax(np.array(mat_grid[sample])))
    if max(maxlist) > 3: # (2)
        # minimum = int(min(minlist)) #- 2
        # maximum = int(max(maxlist)) #+ 2
        minimum = min(minlist) #- 0.025
        maximum = max(maxlist) #+ 0.025
    else:
        minimum = min(minlist) #- 0.025
        maximum = max(maxlist) #+ 0.025
    
    # !!! I have hashed the adjustments out, not sure they are required...
    # (1) this creats a list with the max/min of each material, ignoring any 
    #     numpy NaNs
    # (2) this is to accomodate for the Fano parameter, integer values give
    #     a nicer colorbar scale. The additions ensure the top range of the  
    #     cbar is not used, I think this makes things clearer. 
    
    
    #-------------------------------plot data---------------------------------
    
    #with placeholder = temp_rcParams({'text.usetex': True}):#, 'font.family': 'sans-serif', 'font.sans-serif': 'Helvetica'}):
    
    fig = plt.figure(figsize=(width,height)) 
    

    gs = GridSpec(nrows, ncols, figure=fig,
                  width_ratios=wr,
                  height_ratios=hr,
                  wspace=wspace, hspace=hspace, 
                  top=top, bottom=bottom, left=left, right=right)

    
    arts = list()
    axs = list()
    
    counter = -1
    for nrow in range(rowcount):
        for ncol in range(colcount):
            counter += 1
            data_name = sample_names[counter] # links sample name to gridspec coord
            
            masked_array = np.ma.array(mat_grid[data_name],
                                       mask=np.isnan(mat_grid[data_name])) # (1)
            
            cmap = copy.copy(mpl.cm.inferno) # (2) inferno #gnuplot
            cmap.set_bad('grey', 1.) # (3)
            
            #===============highlight specific point for reference============
            # xfind = 80 # for first map
            # yfind = 25
            # if counter == 1: # for second map
            #     xfind = 90
            #     yfind = 30

            # mat_grid[data_name][xfind][yfind] = np.nan
            # masked_array = np.ma.array(mat_grid[data_name],
            #                             mask=np.isnan(mat_grid[data_name]))
            # cmap.set_bad('c', 1.)
            #=================================================================
            
            # (1) create a masked array, that overlays a boolean array that  
            #     gives True if condition np.isnan is met
            # (2) create a copy of a colour scale
            # (3) if colour scale encounteres True/"1" then set to megenta
            
            #maximum = 12.909584043004712
            #minimum = 2.990106485441589e-32
            
            ax = plt.subplot(gs[nrow,ncol]) # (1)
            art = ax.imshow(masked_array, vmin = minimum, vmax = maximum,
                            cmap=cmap, interpolation=inter, aspect = 'equal')  # (2)
            if log == True:
                if minimum == 0:
                    minimum = 1e-20 # stops log0 error
                # if maximum == np.inf:
                #     maximum = 1e20
                art = ax.imshow(masked_array, cmap=cmap, interpolation=inter, 
                                aspect = 'equal',  #or auto
                                norm=mpl.colors.LogNorm(vmin=minimum,vmax= maximum))
        
                
            axs.append(ax) # append ax and art to their repective lists,
            arts.append(art) # these will have the same index
            
            # (1) save gridspec location in list 'ax' 
            # (2) save map image in list 'art'
            
            
    # add colourbar to last gridspec axis        
    #cb = plt.colorbar(arts[0], cax=plt.subplot(gs[:,-1]))
    if orientation == 'v':
        cb = plt.colorbar(arts[0], cax=plt.subplot(gs[-1,:]), 
                          orientation='horizontal')
    elif orientation == 'h':
        cb = plt.colorbar(arts[0], cax=plt.subplot(gs[:,-1]))
    #cb = plt.colorbar(arts[0], location = 'bottom', orientation = 'horizontal', fraction=0.01)
    cb.ax.tick_params(labelsize=textsize, zorder=10)
    cb.set_label(axis_label, fontsize=textsize)
    #cb.ax.yaxis.set_tick_params(color='k')
    #cb.ax.set_yticklabels([])
    #cb.ax.set_xticklabels(cb.ax.get_xticklabels(), rotation=45)

    # #=========Change the number of ticks on the logscale colorbar==========
    # num_ticks = 4  # Specify the desired number of ticks
    # # Calculate the tick values and labels
    # #ticks = np.logspace(np.log10(data.min()), np.log10(data.max()), num_ticks)
    # ticks = np.array([0.00001, 0.001, 0.1, 10])
    # tick_labels = ['$10^-5$', '$10^-3$', '$10^-1$', '$10^+1$']#['{:.0e}'.format(t) for t in ticks]  # Format tick labels as desired
    # #tick_labels = [tl.replace('1e', '10^') for tl in tick_labels]
    # cb.locator = FixedLocator(ticks)  # Set the locator to use the specified tick values
    # cb.formatter = FixedFormatter(tick_labels)  # Set the formatter to use the specified tick labels
    # # Update the colorbar
    # cb.update_ticks()
    # #======================================================================      
    
    cnt = 0
    for ax in axs:
        ax.set_xticks([])
        ax.set_yticks([])
    # # ========= add scalebar to all maps =========
    #     # scalebar start position
    #     xpos = ax.get_xlim()[0]+ 0.5 + (1-bar_length[cnt]) # might need to change
    #     ypos = ax.get_ylim()[-1] + 0.5 + 1 + 0.8*bar_height[cnt]
    #     pos = (xpos, ypos)
    #     ax.add_patch(Rectangle(pos, bar_length[cnt], bar_height[cnt],
    #                             color='k', lw = 0, clip_on=False,
    #                             transform=ax.transAxes))
    #     cnt += 1
    #     #ax.set_title(variable, pad = 10.325)
    # # ===========================================
    
    #======= add one scalebar to given map =======
    ax = axs[0]#axs[2]
    xpos = ax.get_xlim()[0]+ 0.5 #+ (1-bar_length[cnt]) # might need to change
    ypos = ax.get_ylim()[-1] + 0.5 + 1 + 0.5*bar_height[cnt] + 0.03
    pos = (xpos, ypos)
    ax.add_patch(Rectangle(pos, bar_length[cnt], bar_height[cnt]+0.01,
                                color='k', lw = 0, clip_on=False,
                                transform=ax.transAxes))
    #=============================================

    # fig.text(xpos-0.2,1.03, str(scalebar_microns)+r'$\mu$m -', 
    #           transform=ax.transAxes) #0.57
    
    # add label to images, transfrom makes numbers represent grid type coords
    #alphabet = ['a','b','c','d','e','f']
    #alphabet = ['i','ii','iii','iv','v','vi']
    #alphabet = [title]
    #alphabet = ['Raw Material', 'Charged', 'Quenched']
    #alphabet = ['Raw', 'Charged', 'Quenched']
    if variable == 'red_Chi2':
        variable = 'red Chi2'
    #alphabet = [sample_names[0], sample_names[1]]#,
               # sample_names[2]]#, sample_names[3]]#, sample_names[4]]
    # for j in range(len(sample_names)):
    #     fig.text(0,1.05, r' ' + alphabet[j], 
    #               transform=axs[j].transAxes)
    #     fig.text(0,1.0325, 'chi: ' + str(thresh_Chi) + ', sigma: ' + str(thresh_sig), 
    #              transform=axs[j].transAxes)
        # fig.text(xpos-0.2,1.0325, str(scalebar_microns)+r'$\mu$m -', 
        #           transform=axs[j].transAxes) #0.57
        
    # cb.ax.yaxis.set_major_locator(ticker.MultipleLocator(50))
    # cb.ax.yaxis.set_minor_locator(ticker.MultipleLocator(25))  
      
    #gs.tight_layout(fig,pad=1, h_pad=None, w_pad=None,)
    if savefigures == True:
        fig.savefig('C:/Users/'+user+'/OneDrive - University College London' + \
                    '/PhD/Python/Data/Figures/' + '_'.join(sample_names) + \
                        '_' + variable + '_cmap.png',
                        transparent = False, dpi=600)#, bbox_inches='tight')
        fig.savefig('C:/Users/'+user+'/OneDrive - University College London' + \
                    '/PhD/Python/Data/Figures/' + '_'.join(sample_names) + \
                        '_' + variable + '_cmap.svg',
                        transparent = True, dpi=600)#, bbox_inches='tight')

    
#-----------------------------end process timer-------------------------------
# tempdp = np.array(tempdp)
# tempdp = tempdp[np.where(np.isnan(tempdp) == False)]

# End process timer
end_time = time.perf_counter()
print("\nScript runtime: %.2f \bs" % (end_time - start_time))

# last runtime = 10.21s

#---------------------------------script end----------------------------------