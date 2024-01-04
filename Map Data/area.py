# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 11:19:58 2021

@author: hlancaster~
"""

from package.Raman import * # imports everything in the Raman package 
                            # this includes various modules and functions 

# lor = lambda x:  I0b*(Gammab/2)**2/((x-cenb)**2+(Gammab/2)**2) + (m * x + c)
# bas = lambda x: m*x + c
        
peaktype_a = 'Lorentzian'    
peaktype_b = 'BWF'

peaka = 1
peakb = 2

idx_A = 'Area ' + str(peaka) + 'ovr' + str(peakb)
idx_I = 'Intensity ' + str(peaka) + 'ovr' + str(peakb)

# add a third peak to sum (divides peak A by peak B+C)
sumpeaks = False
if sumpeaks == True:
    peaktype_c = 'Lorentzian'
    peakc = 2
    idx_A = 'Area ' + str(peaka) + 'ovr' + str(peakb) + str(peakc)
 
import scipy.integrate as integrate
import numpy as np
import pandas as pd
            
for sample in sample_names: 
    
    # if sample == sample_names[1]:
    #     peaka = 4
    #     peakb = 6
    # else:
    #     peaka = 1
    #     peakb = 3
    
    area_df = pd.DataFrame(index=[idx_A])
    int_df = pd.DataFrame(index=[idx_I])
    
    coord = list(clean_data[sample])[0]
    #print(coord)
    xs = clean_data[sample][coord][2,:] # get clean data
    intmin = min(xs)
    intmax = max(xs)
    xmin, xmax, ymin, ymax, xstep, ystep = mapinfo(raw_data, sample)
    # print(xmin)
    # print(xmax)
    # print(xstep)
    
    for x in np.arange(xmin, xmax+xstep, xstep):
        x = np.round(x, decimals = 2)
        for y in np.arange(ymin, ymax+ystep, ystep):  
            y = np.round(y, decimals = 2)
            
            # x = -6#18
            # y = -6#8
            fd = fit_data[sample][x,y]
            
            n = peaka
            i = (n-1)*3
            I0a = fd[i]
            cena = fd[i+1]
            Gammaa = fd[i+2]
            
            n = peakb
            i = (n-1)*3
            I0b = fd[i]
            cenb = fd[i+1]
            Gammab = fd[i+2]
            
            if sumpeaks == True:
                n = peakc
                i = (n-1)*3
                I0c = fd[i]
                cenc = fd[i+1]
                Gammac = fd[i+2]                
            
            m = fd.loc['m']
            c = fd.loc['c']
            
            if np.isnan(cena) == False:
            
                if peaktype_a == 'Lorentzian':
                    inta_temp, abserra = integrate.quad(Lorentzian, intmin, intmax,
                                                        args=(I0a, cena, Gammaa,
                                                              m, c))
                elif peaktype_a == 'BWF':
                    q = fd.loc['q']
                    inta_temp, abserra = integrate.quad(BWF, intmin, intmax,
                                                        args=(I0a, cena, Gammaa,
                                                              q, m, c))            
                
                if peaktype_b == 'Lorentzian':
                    intb_temp, abserrb = integrate.quad(Lorentzian, intmin, intmax,
                                                        args=(I0b, cenb, Gammab,
                                                              m, c))
                elif peaktype_b == 'BWF':
                    q = fd.loc['q']
                    intb_temp, abserrb = integrate.quad(BWF, intmin, intmax,
                                                        args=(I0b, cenb, Gammab,
                                                              q, m, c))                        
                
                intbaseline, abserrb = integrate.quad(Baseline, intmin, 
                                                      intmax, args=(m, c))
                
                inta = inta_temp - intbaseline
                intb = intb_temp - intbaseline
                
                if sumpeaks == True:
                    if peaktype_c == 'Lorentzian':
                        intc_temp, abserrc = integrate.quad(Lorentzian, intmin, intmax,
                                                            args=(I0c, cenc, Gammac,
                                                                  m, c))
                    elif peaktype_c == 'BWF':
                        q = fd.loc['q']
                        intc_temp, abserrc = integrate.quad(BWF, intmin, intmax,
                                                            args=(I0c, cenc, Gammac,
                                                                  q, m, c))     
                    intc = intc_temp - intbaseline
                
                if intb != 0:
                    Aratio = inta/intb
                    if inta == 0: # if numerator zero, set as v small value
                        Aratio = 1e-18 # instead (avoids log0 issues later)
                    if sumpeaks == True:
                        Aratio = inta/(intb+intc)
                        #Aratio = (inta+intb)/intc
                else:
                    Aratio = np.nan
                    print(x,y)

            
            else:
                Aratio = np.nan
                print('t1: ' + str(x) + ' ' + str(y))
            
            
            
            
            area_df[x,y] = Aratio
            # Iratio = I0a/I0b
            # int_df[x,y] = Iratio
            if sumpeaks != True:
                Iratio = I0a/I0b
                int_df[x,y] = Iratio
            
    if (fit_data[sample].index == idx_A).any() == True:
        fit_data[sample] = fit_data[sample].drop(idx_A,axis=0) # for reruns only
        fit_data[sample] = fit_data[sample].drop(idx_I,axis=0)
    
    fit_data[sample] = pd.concat([fit_data[sample], area_df])
    #fit_data[sample] = pd.concat([fit_data[sample], int_df]) # is not the sum peaks value !!!
    if sumpeaks != True:
        fit_data[sample] = pd.concat([fit_data[sample], int_df])
    

