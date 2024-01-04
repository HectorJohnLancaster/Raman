#--------------------------------file notes-----------------------------------

# Creates a combined plot containing all map height histograms and their 
# normal curves. 

# To do:
    # Sort out print info at end, currently it's info overload

#------------------------------import modules---------------------------------

import matplotlib.pyplot as plt
import numpy as np
import time
import scipy.stats as stats
import matplotlib.ticker as ticker

 
#--------------------------------user inputs----------------------------------


median_filter = False
savefigures = True
#fit_var = ['cen2']
# Choose if figure is half or full page width
full_width = 6.75
half_width = 3.375
cmap_height = height
width = 6.5#3.25#3.12#full_width #4.25#cmap_height # change this part
#height = 0.32 * width # can set aspect ratio here
height = 2#1.35#1.85#1.4#1.3#

# Choose the order in which the plots appear, higher number = appear front. 
# Hash out any extra references. 
zorder = dict()
order = [0.3,0.1,0.6,0.9,0.8,0.63,4,5,6,10,8,8,8,8,8,8]
o = 0
for sample in sample_names:
    zorder[sample] = order[o] 
    o += 1
    
# fit_var = ['cen1', 'Gamma1', 'cen2', 'Gamma2', 'cen3', 'Gamma3', 'q',
#         'Area 2ovr3', 'Area 1ovr3']
#fit_var = ['cen3']
#annotation = ['Raw', 'KC$_{10}$', 'Discharged']#'NaC$_{10}$']
annotation = [sample_names[0]]#, sample_names[1], sample_names[2]]#, sample_names[3]]#, sample_names[4]]
#----------------------------start process timer------------------------------

start_time = time.perf_counter()
    
#-----------------------------process fit data--------------------------------

dataset = dict() # initialise variable
mean_data = pd.DataFrame(index=fit_var, columns=sample_names)

#sample_names = [sample_names[0],sample_names[1]]#,sample_names[-1]]

for variable in fit_var:
    for sample in sample_names:
    
        loc = 'best'
        
        temp = mat_grid[sample].values # gets values from heatmap
        t = list() # temporary list
        for i in temp: # have to dig into the heatmap grid data structure
            for j in i: # to get all values in one list
                t.append(j)
        temp = np.array(t) # make array for operations
        temp = temp[np.logical_not(np.isnan(temp))] # get rid of any nans

        if variable == 'q':
            dataset[sample] = temp#abs(1/temp) # (5)
            axis_label = 'abs. Fano factor $| 1/q |$'
        elif variable[:5] == 'Gamma':
            dataset[sample] = temp
            axis_label = 'FWHM G (cm$^{-1}$)'
        elif variable[:3] == 'Int':
            dataset[sample] = temp
            p1 = str(variable[4])
            p2 = str(variable[8])
            axis_label = 'I$_'+p1+'$/I$_{'+p2+'}$' 
            axis_label = 'I$_D$/I$_G$'
        elif variable[:4] == 'Area':
            dataset[sample] = temp
            p1 = str(variable[5])
            p2 = str(variable[9])
            if sumpeaks == True:
                p2 = str(variable[9])+'+'+str(variable[10])
            axis_label = 'A$_'+p1+'$/A$_{'+p2+'}$'
            axis_label = 'A$_D$/A$_G$'
            #loc = 'upper right'
        elif variable[:2] == 'I0':
            dataset[sample] = temp
            axis_label = 'Intensity' 
        # elif variable[:4] == 'cen1':
        #     dataset[sample] = temp
        #     axis_label = 'Intensity'
        #     loc = 'upper right'
        else:
            dataset[sample] = temp
            axis_label = 'Raman shift (cm$^{-1}$)'
            axis_label = 'PosG (cm$^{-1}$)'
            
    #axis_label = 'I$_{LFP}$/I$_G$'
    # (1) select material dataframe
    # (2) create a temporary file w/ chosen variable data as a np array
    # (3) delete first item in temp, this is a string w/ the variable name
    # (4) set data type for manipulation
    # (5) this is the fano factor, 1/q, take the absolute value

    
    x0, sigma = stats.norm.fit(dataset[sample])  # stats


    
    
    #---------------------------define boundaries---------------------------------
    
    maxlist = list()
    minlist = list()
    
    # find max and min over all modified datasets to set the histogram range
    for sample in sample_names:    
        minlist.append(np.amin(dataset[sample]))
        maxlist.append(np.amax(dataset[sample]))
    
    if variable == 'q':
        minimum = round((min(minlist)), 2)
        maximum = round((max(maxlist)), 2)
        bins = int((maximum - minimum)*100 + 1)  
    elif variable[:3] == 'Int':
        minimum = round((min(minlist)), 1)# - 0.1
        maximum = round((max(maxlist)), 1)# + 0.1 
        bins = int((maximum - minimum)*20 + 1)
    elif variable[:4] == 'Area':
        minimum = round((min(minlist)), 2)#- 0.1#- 0.04
        maximum = round((max(maxlist)), 2)# - 1#+ 0.04
        bins = int((maximum - minimum)*30 + 1)
    elif variable[:5] == 'Gamma':
        minimum = round(int(min(minlist)), -1) - 5
        maximum = round(int(max(maxlist)), -1) + 5
        bins = (maximum - minimum)//2 + 1
    elif variable[:3] == 'cen':
        minimum = round(min(minlist), 1)#- 2
        maximum = round(max(maxlist), 1) #+ 2
        bins = int((maximum - minimum)/2) + 1
    elif variable[:2] == 'I0':
        minimum = round(int(min(minlist)), -1)
        maximum = round(int(max(maxlist)), -1) + 10   
        bins = int((maximum - minimum)/5 + 1)


    #bins = 25 
    if bins == 1:
        minimum = round((min(minlist)), 2)
        maximum = round((max(maxlist)), 2)
        bins = int((maximum - minimum)*100 + 1)        
           
    if bins < 10:
        print('y')
        bins = bins*4 - 3

    bins = 70#29
    binps = np.linspace(minimum,maximum,bins)
    
    
    #---------------------------------plotting------------------------------------
    
    # plot defaults
    #colors = ['cornflowerblue', 'plum', 'mediumseagreen', 'lightcoral']
    #colors = ['darkorchid','darkorange']
    #colors = ['mediumaquamarine', 'lightcoral']
    #colors = ['tab:olive', 'tab:red']
    #colors = ['red', 'blue', 'green', 'darkorange', 'c', 'gold','m','orange','pink']
    colors = ['gold', 'purple', 'lightcoral']
    #colors.reverse()
    #colors = ['darkorange']
    #colors = ['gold', 'purple','forestgreen', 'darkorange','lightcoral'] # viridis
    cmap = mpl.cm.get_cmap('inferno')
    colors = [cmap(0.9), cmap(0.1),cmap(0.55),cmap(0.2),cmap(0.65),cmap(0.6),cmap(0.7),cmap(0.8),cmap(0.9)]
    colors = [cmap(0.1), cmap(0.8),cmap(0.2),cmap(0.9),cmap(0.3),cmap(0.6),cmap(0.7),cmap(0.8),cmap(0.9)]
    #colors = [cmap(0.85), cmap(0.1),cmap(0.75),cmap(0.4),cmap(0.5),cmap(0.6),cmap(0.7),cmap(0.8),cmap(0.9)]
    colors = [cmap(0.1), cmap(0.3),cmap(0.5),cmap(0.7),cmap(0.3),cmap(0.6),cmap(0.7),cmap(0.8),cmap(0.9),
              cmap(0.1), cmap(0.3),cmap(0.5),cmap(0.7),cmap(0.3),cmap(0.6),cmap(0.7),cmap(0.8),cmap(0.9)]
    c2 = '#56B4E9'
    c3 = '#009E73'
    c4 = '#F0E442'
    c5 = '#CC79A7'
    c1 = '#E69F00'
    colors = [c2,c1,c3,c4,c5]
    lines = ['-', '--', '-.', ':'] # list of linestyles
    alphas = [0.7,0.8,0.9,0.7, 0.6, 0.6, 0.6, 0.6,0.6,0.6,0.9,1,1,1,1,1,1,1,1,1] # list of transparancies
    
    # initialise parameters
    gaussian = dict()
    histogram = dict()
    
    counter = -1
    
    #with temp_rcParams({'text.usetex': True}):#, 'font.family': 'sans-serif', 'font.sans-serif': 'Helvetica'}):
    # initialise the figure, it's dimensions and axis parameters 'ax'
    fig, ax  = plt.subplots(figsize=(width,height))
    
    x0 = dict()
    sigma = dict()
    labels = list()
    s_err = dict()
    for sample in sample_names:
        
        labels.append(sample)
        
        #print(np.mean(dataset[sample]))
        x0[sample], sigma[sample] = stats.norm.fit(dataset[sample]) # get mean and standev   
        # create an evenly spaced array of 100 point between dataset range.
        xs = np.linspace(minimum, maximum, 100)
        gaussian[sample] = stats.norm.pdf(xs,x0[sample],sigma[sample]) # create Gaussian fit data
        s_err[sample] = sigma[sample]/np.sqrt(len(dataset[sample])) # standard error of the mean
        
        
        counter += 1 # see 'histogram notes' below for more details
        hatches = ['\\\\\\','','','','','','','','\\\\\\','\\\\\\']
        hatches = ['','','','','\\\\\\','','\\\\\\','\\\\\\',
                   '','','','','','','','','','','','','','']
        #labels = ['Raw Material', 'Charged', 'Discharged']
        #labels = ['Raw', 'KC$_{10}$', 'Quenched']
        #labels = ['Raw', 'KC$_{24}$', 'KC$_{10}$', 'KC$_8$']
        #labels = ['514.5 nm', '785 nm']
        #labels = ['ALEEES', 'Sample02']#, '10', '8']
        labels = copy.deepcopy(sample_names)
        
        histogram[sample] = plt.hist(dataset[sample], bins=binps,
                                rwidth=0.95, density=True,
                                alpha=alphas[counter], color=colors[counter],
                                zorder = zorder[sample], 
                                hatch = hatches[counter], edgecolor='w',
                                label = labels[counter])

        # plt.plot(xs, gaussian[sample], linestyle=lines[counter], zorder = 10, 
        #           label = labels[counter] + " Gaussian", color = "k")
        #plt.legend(loc=loc, frameon=False) #, fontsize=7)
        plt.xlabel(axis_label, labelpad=1, color = 'k')
        plt.ylabel('Normalised probability\n density (arb. units)')#, rotation=-90, labelpad=13)
        #plt.ylabel('temp', color='k')#, rotation=-90, labelpad=13)
        ax.yaxis.set_label_position("right")
        
        mean_data.loc[variable, sample] = x0[sample]
        # ax.text(x0[sample], max(histogram[sample][0]), annotation[counter],
        #         rotation = 90, ha = 'center', transform=ax.transData)
        
        # if variable == 'q':
        #     ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
        #     ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))
        #     ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))    
        if variable[:2] == 'I0':
            ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(20))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05)) 
        if variable[0:3] == 'cen':
            ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1))   
        if variable[0:5] == 'Gamma':
            ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05)) 
        if variable[:3] == 'Int':
            ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))   
        if variable[:4] == 'Area':
            ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))  
        if variable[:2] == 'I0':
            ax.xaxis.set_major_locator(ticker.MultipleLocator(1000))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(200))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))     

        
            
    # fig.text(0,1.02, 'mean: %.2f' % (x0[sample]) + ' stdv: %.2f' % (sigma[sample]), 
    #           transform=ax.transAxes, color='r') #0.57
        
    #plt.title(variable)
    #2.98
    #ax.set_xlim([1592,1620])
    #plt.yticks(rotation=-90)
    #plt.xticks(rotation=-30)
    
    # ticks = list(ax.get_xticks())
    # ticks = [int(item) for item in ticks]
    # ax.set_xticklabels(ticks)
    #plt.xticks([]) 
    plt.yticks([])
    ax.ticklabel_format(useOffset=False)
    
    plt.xticks(zorder=10)
    
    # histogram notes:
        # density=True sets to a normalised probability density
        # rwidth sets the column width               
        
    # save figure
    fig.tight_layout(pad=0.5)
    
    if savefigures == True:
        #print('t')
        fig.savefig('C:/Users/'+user+'/OneDrive - University College London' + \
                    '/PhD/Python/Data/Figures/' + '_'.join(sample_names) + '_' + \
                    variable + '_hist.png', transparent = True, dpi=600)
        fig.savefig('C:/Users/'+user+'/OneDrive - University College London' + \
                    '/PhD/Python/Data/Figures/' + '_'.join(sample_names) + '_' + \
                    variable + '_hist.svg', transparent = True, dpi=600)
        
    
    #---------------------------------print data----------------------------------
    
    for sample in sample_names:
        print("\n" + sample + ": mean = (%.5f +/- %0.5f) cm-1 and 2sd = %.5f" % \
              (x0[sample], s_err[sample], sigma[sample]*2))
    
    
    #-----------------------------statistical tests-------------------------------
    
    # All code in this section modified from machinelearningmaster.com
    # "A Gentle Introduction to Normality Tests in Python"
    
    for sample in sample_names:
        data = dataset[sample]
    
        #-----Shapiro-Wilk Test-----
        from scipy.stats import shapiro
        stat, p = shapiro(data)
        print('\n---' + sample + '---\nShapiro-Wilk Statistic: %.3f, p=%.3f' % (stat, p))
        # interpret
        alpha = 0.05
        if p > alpha:
         	print('Looks Gaussian (fail to reject H0)')
        else:
         	print('Does not look Gaussian (reject H0)')
            
        #-----D'Agostino and Pearson Test----- 
        from scipy.stats import normaltest
        stat, p = normaltest(data)
        print('\nDAgostino and Pearson Statistic: %.3f, p=%.3f' % (stat, p))
        # interpret
        alpha = 0.05
        if p > alpha:
         	print('Looks Gaussian (fail to reject H0)')
        else:
         	print('Does not look Gaussian (reject H0)')
            
        #-----Anderson-Darling Test-----
        from scipy.stats import anderson
        result = anderson(data)
        print('\nAnderson-Darling Statistic: %.3f' % result.statistic)
        p = 0
        for i in range(len(result.critical_values)):
         	sl, cv = result.significance_level[i], result.critical_values[i]
         	if result.statistic < result.critical_values[i]:
                 print('%.3f: %.3f, looks Gaussian (fail to reject H0))' % (sl, cv))
         	else:
                 print('%.3f: %.3f, does not look Gaussian (reject H0)' % (sl, cv))

# nrows=2
# ncols=4
# fig = plt.figure(figsize=(width,height)) 
# gs = GridSpec(nrows, ncols, figure=fig,
#               left=0.025, right=0.925, top=0.9, bottom=0.05,
#               wspace=0.51)    
    
# fig, axs = plt.subplots(nrows=4, ncols=2, constrained_layout=True, figsize=(8,8))
# counter = 0
# condition = True
# titles = ['cen D', 'FWHM D', 'cen G', 'FWHM G', 'abs. Fano', 'I$_D$/I$_G$',
#           'A$_D$/A$_G$']

# while(condition):
#     for ax in axs.flat:
#         var = fit_var[counter]
#         title = titles[counter]
#         ytest = np.array(mean_data.loc[var])
#         xtest = ['Raw', 'Charge', 'Q1h', 'Q2w', 'Q4w']
#         ax.plot(xtest,ytest, marker = 'D', color = 'k', linestyle=':')
#         ax.set_title(title)
#         counter += 1
#         if counter == 7:
#             condition = False
#             break
# fig.savefig('C:/Users/hlancaster/OneDrive - University College London' + \
#             '/PhD/Python/Data/Figures/change.tif', dpi=1200)
                
#-----------------------------end process timer-------------------------------

# End process timer
end_time = time.perf_counter()
print("\nScript runtime: %.2f \bs" % (end_time - start_time))
# last runtime = 3.31s

#---------------------------------script end----------------------------------