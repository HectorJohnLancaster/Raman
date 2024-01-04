#--------------------------------file notes-----------------------------------

# Applies a Lorentzian fit to the D peak and Breit-Wigner-Fano fit to the G 
# peak. 

# To do:
#       4) Quantatative comparrison of fits

#--------------------------------user inputs----------------------------------

# Choose which data to chuck, below this threshold of reduced Chi2 (thresh_Chi)
# and fit_var standard deviations (thresh_sig) is kept.
pickle_data = False
thresh_Chi = 20 #20 default
thresh_sig = 6 #6 default

fit_type_list = ['2Lor', '3Lor', '4Lor', '5Lor', '9Lor',
                 'LorBWF', '2LorBWF', '3LorBWF', '4LorBWF',
                 '2Voigt', '3Voigt', '4Voigt', '5Voigt']
material_list = ['Active Carbon', 'MWCNT', 'Carbon Electrode 1',
                 'Carbon Electrode 2', 'Carbon Electrode 3', 
                 'Carbon Electrode 4', 'Carbon Electrode 5', 'MnOx', 
                 'Alex Diamond', 'Black Phosphorus', 'CD Grain 7 3L',
                 'CD Grain 7 4L', 'CD Grain 6 2L', 'CD Grain 6 3L',
                 'CD Grain 6 4L', 'CD Grain 6 5L', 'CD Grain 6 2V',
                 'CD Grain 6 3V', 'CD Grain 6 4V', 'CD Grain 6 5V',
                 'CD Grain 8 2L', 'CD Grain 8 3L', 'CD Grain 8 4L']

fit_type_list = ['1Lor1BWF1Lor']#['2Lor1BWF', '2Lor1BWF', '2Lor1BWF']#, '2Lor1BWF', '2Lor1BWF']
#fit_type_list = list()
# for name in sample_names:
#     fit_type_list.append('2Lor1BWF')
#     # # if name[:2] == 'CB':
#     # #     fit_type_list.append('3LorBWF')
#     # if name[:2] == 'CB' or name[3:5] =='CB':
#     #     fit_type_list.append('2Lor1BWF1Lor')
#     # else:
#     #     fit_type_list.append('2Lor1BWF')
        
material_type_list = ['SC Graphite']#['Active Carbon', 'Active Carbon', 'Active Carbon']#['Carbon Black 4p', 'Carbon Black 4p'] #'Active Carbon', 'Active Carbon']#, 'Active Carbon', 'Active Carbon']
# material_type_list = list()
# for name in sample_names:
#     material_type_list.append('Active Carbon')
#     # # if name == 'AC kc8 map':
#     # #     material_type_list.append('AC kc8')
#     # if name[:2] == 'CB' or name[3:5] =='CB':
#     #     material_type_list.append('Carbon Black 4p')
#     # else:
#     #     material_type_list.append('Active Carbon')
        

#-----------------------------import modules----------------------------------

import pickle
import time
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt
import sys
import scipy.stats as stats
import warnings
import initial_guesses # Must be in the working directory
from warnings import simplefilter
import winsound
frequency = 500  # Set Frequency To 2500 Hertz
duration = 1000  # Set Duration To 1000 ms == 1 second
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)


#----------------------------start process timer------------------------------

start_time = time.perf_counter()

#-----------------------------define functions--------------------------------


#_____________________________individual curves_______________________________

# See cauchy distribution function at:
# https://www.itl.nist.gov/div898/handbook/eda/section3/eda3663.htm 
# Here, the  modified version with the max intensity (I0D) replaces the 
# the parameter (1/s*pi).
def Lorentzian(x, I0, cen, Gamma, m, c):
    return I0*(Gamma/2)**2/((x-cen)**2+(Gamma/2)**2) + (m * x + c)

# Function as described in SI of DOI: 10.1103/PhysRevB.84.241404    
def BWF(x, I0, cen, Gamma, q, m, c):
    numerator = (1+(x-cen)/(q*(Gamma/2)))**2
    denominator = 1+((x-cen)/(Gamma/2))**2
    return I0*(numerator/denominator) + (m * x + c)

# Linear function used as baseline
def Baseline(x, m, c):
    return m * x + c

def spectrum_interval(x,y, lb, ub):
    lowerx = np.where(x>lb)
    upperx = np.where(x<ub)
    interx = np.intersect1d(lowerx,upperx)
    x = x[interx]  
    y = y[interx]
    return(x,y)

# Extract information about the map setup
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
    return ltrim[1] - ltrim[0] # return a step size


    
    
    
#--------------------------get & print fitting data---------------------------

plt.rcParams.update({'figure.max_open_warning': 0}) # supresses runtime info


fit_data = dict() # initialise
report = dict()
s_count = 0
sample_num = str(len(sample_names))
cnt = 0
for sample in sample_names: # for each map

    #_____________________________combined function___________________________
    
    fit_type = fit_type_list[cnt]
    material_type = material_type_list[cnt]

    if fit_type == '1Lor':
         def ComFit(x, I01, cen1, Gamma1, m, c):
             Lor1 = Lorentzian(x, I01, cen1, Gamma1, m, c)
             return Lor1
    
         # Choose parameter names 'p_names' for the dataframe, order must match  
         # the order that the parameters are storred upon entry.
         p_names = ['I01', 'cen1', 'Gamma1',
                    'm', 'c', 'Chi2','red_Chi2']

    elif fit_type == '2Lor':
         def ComFit(x, I0D, cenD, GammaD, I0G, cenG, GammaG, m, c):
             Lor1 = Lorentzian(x, I0D, cenD, GammaD, m, c)
             Lor2 = Lorentzian(x, I0G, cenG, GammaG, m, c)
             return Lor1 + Lor2 - (m * x + c)
    
         p_names = ['I01', 'cen1', 'Gamma1', 
                    'I02', 'cen2', 'Gamma2', 
                    'm', 'c', 'Chi2','red_Chi2']
         
         
    elif fit_type == '1Lor1BWF':
         def ComFit(x, I01, cen1, Gamma1, I02, cen2, Gamma2, q, m, c):
             Lor1 = Lorentzian(x, I01, cen1, Gamma1, m, c)
             BWF1 = BWF(x, I02, cen2, Gamma2, q, m, c)
             return Lor1 + BWF1 - (m * x + c) 
    
         p_names = ['I01', 'cen1', 'Gamma1', 
                    'I02', 'cen2', 'Gamma2', 
                    'q', 'm', 'c', 'Chi2','red_Chi2']
        
        
    elif fit_type == '2Lor1BWF':
        def ComFit(x, I01, cen1, Gamma1, 
                   I02, cen2, Gamma2,
                   I03, cen3, Gamma3, 
                   q, m, c):
            Lor1 = Lorentzian(x, I01, cen1, Gamma1, m, c)
            Lor2 = Lorentzian(x, I02, cen2, Gamma2, m, c)
            BWF3 = BWF(x, I03, cen3, Gamma3, q, m, c)
    
            baseline = m * x +c
            return Lor1 + Lor2 + BWF3 - 2*baseline
        
        p_names = ['I01', 'cen1', 'Gamma1', 
                   'I02', 'cen2', 'Gamma2',
                   'I03', 'cen3', 'Gamma3', 
                   'q', 'm', 'c', 'Chi2','red_Chi2']
        
    elif fit_type == '3Lor1BWF':
        def ComFit(x, I01, cen1, Gamma1, 
                    I02, cen2, Gamma2,
                    I03, cen3, Gamma3,
                    I04, cen4, Gamma4,
                    q, m, c):
            Lor1 = Lorentzian(x, I01, cen1, Gamma1, m, c)
            Lor2 = Lorentzian(x, I02, cen2, Gamma2, m, c)
            Lor3 = Lorentzian(x, I03, cen3, Gamma3, m, c)
            BWF1 = BWF(x, I04, cen4, Gamma4, q, m, c)
            baseline = m * x +c
            return Lor1 + Lor2 + Lor3 + BWF1 - 3*baseline
        
        p_names = ['I01', 'cen1', 'Gamma1', 
                    'I02', 'cen2', 'Gamma2',
                    'I03', 'cen3', 'Gamma3',
                    'I04', 'cen4', 'Gamma4',
                    'q', 'm', 'c', 'Chi2','red_Chi2']

    elif fit_type == '1Lor1BWF1Lor':
        def ComFit(x, I01, cen1, Gamma1, 
                    I02, cen2, Gamma2,
                    I03, cen3, Gamma3,
                    q, m, c):
            Lor1 = Lorentzian(x, I01, cen1, Gamma1, m, c)
            BWF1 = BWF(x, I02, cen2, Gamma2, q, m, c)
            Lor2 = Lorentzian(x, I03, cen3, Gamma3, m, c)
            baseline = m * x +c
            return Lor1 + Lor2 + BWF1 - 2*baseline
        
        p_names = ['I01', 'cen1', 'Gamma1', 
                   'I02', 'cen2', 'Gamma2',
                   'I03', 'cen3', 'Gamma3',
                   'q', 'm', 'c', 'Chi2','red_Chi2']
        
    elif fit_type == '2Lor1BWF1Lor':
        def ComFit(x, I01, cen1, Gamma1, 
                    I02, cen2, Gamma2,
                    I03, cen3, Gamma3,
                    I04, cen4, Gamma4,
                    q, m, c):
            Lor1 = Lorentzian(x, I01, cen1, Gamma1, m, c)
            Lor2 = Lorentzian(x, I02, cen2, Gamma2, m, c)
            BWF1 = BWF(x, I03, cen3, Gamma3, q, m, c)
            Lor3 = Lorentzian(x, I04, cen4, Gamma4, m, c)
            baseline = m * x +c
            return Lor1 + Lor2 + Lor3 + BWF1 - 3*baseline
        
        p_names = ['I01', 'cen1', 'Gamma1', 
                   'I02', 'cen2', 'Gamma2',
                   'I03', 'cen3', 'Gamma3',
                   'I04', 'cen4', 'Gamma4',
                   'q', 'm', 'c', 'Chi2','red_Chi2']
     
    elif fit_type == '3Lor':
        def ComFit(x, I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, m, c):
            Lor1 = Lorentzian(x, I01, cen1, Gamma1, m, c)
            Lor2 = Lorentzian(x, I02, cen2, Gamma2, m, c)
            Lor3 = Lorentzian(x, I03, cen3, Gamma3, m, c)
            return Lor1 + Lor2 + Lor3 - 2*(m * x + c)
        
        p_names = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                   'I03', 'cen3', 'Gamma3', 
                   'm', 'c', 'Chi2','red_Chi2']  
        
        
    elif fit_type == '4Lor':
        def ComFit(x, I01, cen1, Gamma1, I02, cen2, Gamma2,
                      I03, cen3, Gamma3, I04, cen4, Gamma4,
                      m, c):
            Lor1 = Lorentzian(x, I01, cen1, Gamma1, m, c)
            Lor2 = Lorentzian(x, I02, cen2, Gamma2, m, c)
            Lor3 = Lorentzian(x, I03, cen3, Gamma3, m, c)
            Lor4 = Lorentzian(x, I04, cen4, Gamma4, m, c)
            baseline = m * x +c
            return Lor1 + Lor2 + Lor3 + Lor4 - 3*baseline
        
        p_names = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                   'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
                   'm', 'c', 'Chi2','red_Chi2']
    
    elif fit_type == '5Lor':
        def ComFit(x, I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, m, c):
            Lor1 = Lorentzian(x, I01, cen1, Gamma1, m, c)
            Lor2 = Lorentzian(x, I02, cen2, Gamma2, m, c)
            Lor3 = Lorentzian(x, I03, cen3, Gamma3, m, c)
            Lor4 = Lorentzian(x, I04, cen4, Gamma4, m, c)
            Lor5 = Lorentzian(x, I05, cen5, Gamma5, m, c)
            return Lor1 + Lor2 + Lor3 + Lor4 + Lor5 - 4*(m * x + c)
                
        p_names = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                   'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
                   'I05', 'cen5', 'Gamma5', 
                   'm', 'c', 'Chi2','red_Chi2'] 
        
    
    elif fit_type == '6Lor':
        def ComFit(x, I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, I06, cen6, Gamma6, m, c):
            Lor1 = Lorentzian(x, I01, cen1, Gamma1, m, c)
            Lor2 = Lorentzian(x, I02, cen2, Gamma2, m, c)
            Lor3 = Lorentzian(x, I03, cen3, Gamma3, m, c)
            Lor4 = Lorentzian(x, I04, cen4, Gamma4, m, c)
            Lor5 = Lorentzian(x, I05, cen5, Gamma5, m, c)
            Lor6 = Lorentzian(x, I06, cen6, Gamma6, m, c)
            return Lor1 + Lor2 + Lor3 + Lor4 + Lor5 + Lor6 - 5*(m * x + c)
                
        p_names = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                   'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
                   'I05', 'cen5', 'Gamma5', 'I06', 'cen6', 'Gamma6',
                   'm', 'c', 'Chi2','red_Chi2']   
        
        
    elif fit_type == '4Lor1BWF':
        def ComFit(x, I01, cen1, Gamma1, I02, cen2, Gamma2,
                      I03, cen3, Gamma3, I04, cen4, Gamma4,
                      I05, cen5, Gamma5, q, m, c):
            Lor1 = Lorentzian(x, I01, cen1, Gamma1, m, c)
            Lor2 = Lorentzian(x, I02, cen2, Gamma2, m, c)
            Lor3 = Lorentzian(x, I03, cen3, Gamma3, m, c)
            Lor4 = Lorentzian(x, I04, cen4, Gamma4, m, c)
            BWF1 = BWF(x, I05, cen5, Gamma5, q, m, c)
            return Lor1 + Lor2 + Lor3 + Lor4 + BWF1 - 4*(m * x + c)
        
        p_names = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                   'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
                   'I05', 'cen5', 'Gamma5', 'q', 'm', 'c', 'Chi2','red_Chi2']   
        
        
    elif fit_type == '5Lor1BWF':
        def ComFit(x, I01, cen1, Gamma1, I02, cen2, Gamma2,
                      I03, cen3, Gamma3, I04, cen4, Gamma4,
                      I05, cen5, Gamma5, I06, cen6, Gamma6, q, m, c):
            Lor1 = Lorentzian(x, I01, cen1, Gamma1, m, c)
            Lor2 = Lorentzian(x, I02, cen2, Gamma2, m, c)
            Lor3 = Lorentzian(x, I03, cen3, Gamma3, m, c)
            Lor4 = Lorentzian(x, I04, cen4, Gamma4, m, c)
            Lor5 = Lorentzian(x, I05, cen5, Gamma5, m, c)
            BWF1 = BWF(x, I06, cen6, Gamma6, q, m, c)
            return Lor1 + Lor2 + Lor3 + Lor4 + Lor5 + BWF1 - 5*(m * x + c)
        
        p_names = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                   'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
                   'I05', 'cen5', 'Gamma5', 'I06', 'cen6', 'Gamma6', 'q', 'm', 'c', 'Chi2','red_Chi2']   
    
    
    elif fit_type == '7Lor':
        def ComFit(x, I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, I06, cen6, Gamma6,
                   I07, cen7, Gamma7, m, c):
            Lor1 = Lorentzian(x, I01, cen1, Gamma1, m, c)
            Lor2 = Lorentzian(x, I02, cen2, Gamma2, m, c)
            Lor3 = Lorentzian(x, I03, cen3, Gamma3, m, c)
            Lor4 = Lorentzian(x, I04, cen4, Gamma4, m, c)
            Lor5 = Lorentzian(x, I05, cen5, Gamma5, m, c)
            Lor6 = Lorentzian(x, I06, cen6, Gamma6, m, c)
            Lor7 = Lorentzian(x, I07, cen7, Gamma7, m, c)
            return Lor1 + Lor2 + Lor3 + Lor4 + Lor5 + Lor6 + Lor7 - 6*(m * x + c)
                
        p_names = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                   'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
                   'I05', 'cen5', 'Gamma5', 'I06', 'cen6', 'Gamma6',
                   'I07', 'cen7', 'Gamma7', 'm', 'c', 'Chi2','red_Chi2']  
    
    elif fit_type == '8Lor':
        def ComFit(x, I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, I06, cen6, Gamma6,
                   I07, cen7, Gamma7, I08, cen8, Gamma8,
                   m, c):
            Lor1 = Lorentzian(x, I01, cen1, Gamma1, m, c)
            Lor2 = Lorentzian(x, I02, cen2, Gamma2, m, c)
            Lor3 = Lorentzian(x, I03, cen3, Gamma3, m, c)
            Lor4 = Lorentzian(x, I04, cen4, Gamma4, m, c)
            Lor5 = Lorentzian(x, I05, cen5, Gamma5, m, c)
            Lor6 = Lorentzian(x, I06, cen6, Gamma6, m, c)
            Lor7 = Lorentzian(x, I07, cen7, Gamma7, m, c)
            Lor8 = Lorentzian(x, I08, cen8, Gamma8, m, c)
            return Lor1 + Lor2 + Lor3 + Lor4 + Lor5 + Lor6 + Lor7 + Lor8 - \
                7*(m * x + c)
                
        p_names = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                   'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
                   'I05', 'cen5', 'Gamma5', 'I06', 'cen6', 'Gamma6',
                   'I07', 'cen7', 'Gamma7', 'I08', 'cen8', 'Gamma8',
                   'm', 'c', 'Chi2','red_Chi2']  
        
    elif fit_type == '9Lor':
        def ComFit(x, I01, cen1, Gamma1, I02, cen2, Gamma2,
                   I03, cen3, Gamma3, I04, cen4, Gamma4,
                   I05, cen5, Gamma5, I06, cen6, Gamma6,
                   I07, cen7, Gamma7, I08, cen8, Gamma8,
                   I09, cen9, Gamma9, m, c):
            Lor1 = Lorentzian(x, I01, cen1, Gamma1, m, c)
            Lor2 = Lorentzian(x, I02, cen2, Gamma2, m, c)
            Lor3 = Lorentzian(x, I03, cen3, Gamma3, m, c)
            Lor4 = Lorentzian(x, I04, cen4, Gamma4, m, c)
            Lor5 = Lorentzian(x, I05, cen5, Gamma5, m, c)
            Lor6 = Lorentzian(x, I06, cen6, Gamma6, m, c)
            Lor7 = Lorentzian(x, I07, cen7, Gamma7, m, c)
            Lor8 = Lorentzian(x, I08, cen8, Gamma8, m, c)
            Lor9 = Lorentzian(x, I09, cen9, Gamma9, m, c)
            return Lor1 + Lor2 + Lor3 + Lor4 + Lor5 + Lor6 + Lor7 + Lor8 + \
                Lor9 - 8*(m * x + c)
                
        p_names = ['I01', 'cen1', 'Gamma1', 'I02', 'cen2', 'Gamma2',
                   'I03', 'cen3', 'Gamma3', 'I04', 'cen4', 'Gamma4', 
                   'I05', 'cen5', 'Gamma5', 'I06', 'cen6', 'Gamma6',
                   'I07', 'cen7', 'Gamma7', 'I08', 'cen8', 'Gamma8',
                   'I09', 'cen9', 'Gamma9',
                   'm', 'c', 'Chi2','red_Chi2']            
    
    
    elif fit_type == '2Voigt':
        def ComFit(x, amp1, sigma1, I01, cen1, Gamma1,
                   amp2, sigma2, I02, cen2, Gamma2, m, c):
            Voigt1 = Voigt(x, amp1, sigma1, I01, cen1, Gamma1, m, c)
            Voigt2 = Voigt(x, amp2, sigma2, I02, cen2, Gamma2, m, c)
            return Voigt1 + Voigt2 - (m * x + c)        
    
        p_names = ['amp1', 'sigma1', 'I01', 'cen1', 'Gamma1',
                   'amp2', 'sigma2', 'I02', 'cen2', 'Gamma2',
                   'm', 'c', 'Chi2','red_Chi2']
    
    elif fit_type == '3Voigt':
        def ComFit(x, amp1, sigma1, I01, cen1, Gamma1,
                   amp2, sigma2, I02, cen2, Gamma2,
                   amp3, sigma3, I03, cen3, Gamma3, m, c):
            Voigt1 = Voigt(x, amp1, sigma1, I01, cen1, Gamma1, m, c)
            Voigt2 = Voigt(x, amp2, sigma2, I02, cen2, Gamma2, m, c)
            Voigt3 = Voigt(x, amp3, sigma3, I03, cen3, Gamma3, m, c)
            return Voigt1 + Voigt2 + Voigt3 - 2*(m * x + c)        
    
        p_names = ['amp1', 'sigma1', 'I01', 'cen1', 'Gamma1',
                   'amp2', 'sigma2', 'I02', 'cen2', 'Gamma2',
                   'amp3', 'sigma3', 'I03', 'cen3', 'Gamma3', 
                   'm', 'c', 'Chi2','red_Chi2']
        
    elif fit_type == '4Voigt':
        def ComFit(x, amp1, sigma1, I01, cen1, Gamma1,
                   amp2, sigma2, I02, cen2, Gamma2,
                   amp3, sigma3, I03, cen3, Gamma3,
                   amp4, sigma4, I04, cen4, Gamma4, m, c):
            Voigt1 = Voigt(x, amp1, sigma1, I01, cen1, Gamma1, m, c)
            Voigt2 = Voigt(x, amp2, sigma2, I02, cen2, Gamma2, m, c)
            Voigt3 = Voigt(x, amp3, sigma3, I03, cen3, Gamma3, m, c)
            Voigt4 = Voigt(x, amp4, sigma4, I04, cen4, Gamma4, m, c)
            return Voigt1 + Voigt2 + Voigt3 + Voigt4 - 3*(m * x + c)        
    
        p_names = ['amp1', 'sigma1', 'I01', 'cen1', 'Gamma1',
                   'amp2', 'sigma2', 'I02', 'cen2', 'Gamma2',
                   'amp3', 'sigma3', 'I03', 'cen3', 'Gamma3', 
                   'amp4', 'sigma4', 'I04', 'cen4', 'Gamma4', 
                   'm', 'c', 'Chi2','red_Chi2']
        
    elif fit_type == '5Voigt':
        def ComFit(x, amp1, sigma1, I01, cen1, Gamma1,
                   amp2, sigma2, I02, cen2, Gamma2,
                   amp3, sigma3, I03, cen3, Gamma3,
                   amp4, sigma4, I04, cen4, Gamma4,
                   amp5, sigma5, I05, cen5, Gamma5, m, c):
            Voigt1 = Voigt(x, amp1, sigma1, I01, cen1, Gamma1, m, c)
            Voigt2 = Voigt(x, amp2, sigma2, I02, cen2, Gamma2, m, c)
            Voigt3 = Voigt(x, amp3, sigma3, I03, cen3, Gamma3, m, c)
            Voigt4 = Voigt(x, amp4, sigma4, I04, cen4, Gamma4, m, c)
            Voigt5 = Voigt(x, amp5, sigma5, I05, cen5, Gamma5, m, c)
            return Voigt1 + Voigt2 + Voigt3 + Voigt4 + Voigt5 - 4*(m * x + c)        
    
        p_names = ['amp1', 'sigma1', 'I01', 'cen1', 'Gamma1',
                   'amp2', 'sigma2', 'I02', 'cen2', 'Gamma2',
                   'amp3', 'sigma3', 'I03', 'cen3', 'Gamma3', 
                   'amp4', 'sigma4', 'I04', 'cen4', 'Gamma4', 
                   'amp5', 'sigma5', 'I05', 'cen5', 'Gamma5', 
                   'm', 'c', 'Chi2','red_Chi2']
        
    
    #-------------------------------
    temp_count = [0,0]
    xmin, xmax, ymin, ymax, xstep, ystep = mapinfo(sample) # get map info
    fit_df = pd.DataFrame(index=p_names) # initialise dataframe
    # --- percentage counter bits ---
    pc_count = 0
    pc_tot = len(np.arange(xmin, xmax+xstep, xstep))
    s_count += 1
    sys.stdout.write('\n')
    # -------------------------------
                
        
    for x in np.arange(xmin, xmax+xstep, xstep):
        x = np.round(x, decimals = 2)
        # --- percentage counter bits ---
        pc_count += 1
        pc = (pc_count*100)//(pc_tot)
        sys.stdout.write('\rFitting ' + str(s_count) + ' of ' + 
                         sample_num + ': ' + str(pc) + '%')
        # -------------------------------
        for y in np.arange(ymin, ymax+ystep, ystep):
            y = np.round(y, decimals = 2)
            temp_count[1] += 1
            #----------------define spectra variables and error---------------  

            xs = clean_data[sample][x,y][2,:]
            ys = clean_data[sample][x,y][3,:]
                
            if material_type == 'Lasered Diaomond Undoped':
                xs, ys = spectrum_interval(xs, ys, 950, 1850)
                # xs = xs[ys != 0] # delete any zeros
                # ys = ys[ys != 0]
            
            if material_type[:4] == 'MnOx':
                xs, ys = spectrum_interval(xs, ys, 200, 850)
                
            if material_type == 'MWCNT':
                xs, ys = spectrum_interval(xs, ys, 1000, 1800)
                
            if material_type == 'Canyon Diablo':
                xs, ys = spectrum_interval(xs, ys, 1000, 1950)
                
            if material_type == 'PtG20-5c':
                xs, ys = spectrum_interval(xs, ys, 1100, 3200)
                
            if sample[-5:] == 'Anode':
                xs, ys = spectrum_interval(xs, ys, 850, 1900)
            
            #xs, ys = spectrum_interval(xs, ys, 700, 1980)    
            
            ys_err = np.sqrt(abs(ys)) # (1)
            ys_err += 1*10**-20 # (2)
            
            # (1) Uncertainty in intensity values, since photons hitting a 
            #     detector follow a Poisson distribution.
            # (2) Adds a tiny uncertainty to avoid dividing by zero later on
            
            
            #---------------------find BWF/Lorentzian fit---------------------   
            
            # The below is taken from the documentation on scipy.org:
            # scipy.optimize.curve_fit(f, xdata, ydata, p0)
            # uses a non-linear least squares to fit a function, f, to data.
            # f = the model function
            # xdata = the independent variable where the data is measured
            # ydata = the dependent data
            # p0 = initial guess for the parameters of f
            # sigma = uncertainty in ydata
            # absolute_sigma = False, set as makes sense with quoted errors
            # method = lm, since the problem is not constrained, this can be
            #   used and offers the fastest computational time. The other two 
            #   methods also work, but take about twice as long. All methods
            #   give the same reduced Chi2 statistic. 
            
            # returns: 
            #   popt: optimal values for the parameters so that the sum of the 
            #         squared residuals of f(xdata, *popt) - ydata is minimized
            #   pcov: the estimated covarience of popt. The diagonals provide 
            #         the variance of the parameter estimate.
            
            
            if np.isnan(ys).all() != True:
    
                checker = True # if the optimization works, checker stays True
                
                try:
                    second_pass = False
                    guesses, bnds = initial_guesses.main(material_type, xs, ys,
                                                         second_pass, sample)
                    popt_ComFit1, pcov_ComFit1 = curve_fit(ComFit, xs, ys,
                                                         p0=guesses,
                                                         sigma = ys_err,
                                                         absolute_sigma=False,
                                                         bounds = bnds,
                                                         method='trf')
                    second_pass = True
                    guesses, bnds = initial_guesses.main(material_type, xs, ys,
                                                          second_pass, sample)
                    popt_ComFit2, pcov_ComFit2 = curve_fit(ComFit, xs, ys,
                                                          p0=guesses,
                                                          sigma = ys_err,
                                                          absolute_sigma=False,
                                                          bounds = bnds,
                                                          method='trf')               
                
                # if optimisation fails, set values to NaN and checker to False
                except RuntimeError:
                    checker = False             
                    # fig, ax = plt.subplots()
                    # ax.plot(xs, ys, label= 'material: ' + \
                    #         sample + '\ncoords: ' + str((x,y)))
                    # ax.legend()
                
                    if second_pass == True:
                        checker = True
                        print('RuntimeError, optimisation at: ' + str((x,y)) + \
                              ' only valid on first pass')
                        popt_ComFit2 = np.zeros(np.shape(guesses))
                        pcov_ComFit2 = np.zeros((len(guesses),len(guesses)))
                        popt_ComFit2[:] = np.nan              
                        pcov_ComFit2[:] = np.nan
    
                    else:
                        print('optimisation at: ' + str((x,y)) + 'failed')
                        popt_ComFit = np.zeros(np.shape(guesses))
                        pcov_ComFit = np.zeros((len(guesses),len(guesses)))
                        popt_ComFit[:] = np.nan              
                        pcov_ComFit[:] = np.nan
                        popt_ComFit1 = popt_ComFit
                        popt_ComFit2 = popt_ComFit
                        pcov_ComFit1 = pcov_ComFit
                        pcov_ComFit2 = pcov_ComFit
    
                except ValueError:
                        print('ValueError, optimisation at: ' + str((x,y)) + 'failed')
                        popt_ComFit = np.zeros(np.shape(guesses))
                        pcov_ComFit = np.zeros((len(guesses),len(guesses)))
                        popt_ComFit[:] = np.nan              
                        pcov_ComFit[:] = np.nan
                        popt_ComFit1 = popt_ComFit
                        popt_ComFit2 = popt_ComFit  
                        pcov_ComFit1 = pcov_ComFit
                        pcov_ComFit2 = pcov_ComFit                     
                
                #2.5
                if 0.5*(max(ys[:10])-min(ys[:10])) > max(ys)-min(ys):
                        print('Data too noisy at: ' + str((x,y)) + ' set to NaN')
                        popt_ComFit = np.zeros(np.shape(guesses))
                        pcov_ComFit = np.zeros((len(guesses),len(guesses)))
                        popt_ComFit[:] = np.nan              
                        pcov_ComFit[:] = np.nan
                        popt_ComFit1 = popt_ComFit
                        popt_ComFit2 = popt_ComFit  
                        pcov_ComFit1 = pcov_ComFit
                        pcov_ComFit2 = pcov_ComFit     
                        noise += 1

                #-----------------------------stats-------------------------------
                
                residual = ys - (ComFit(xs, *popt_ComFit1)) # (1)
                normres = residual/ys_err # (2)            
                Chi2 = sum(normres**2) # (3)
                K = len(xs) - len(popt_ComFit1) # (4)
                red_Chi2 = Chi2/K # (5)
                
                residual = ys - (ComFit(xs, *popt_ComFit2)) 
                normres = residual/ys_err          
                sChi2 = sum(normres**2)
                K = len(xs) - len(popt_ComFit2)
                sred_Chi2 = sChi2/K 
                
                if sred_Chi2 < red_Chi2:
                    popt_ComFit = popt_ComFit2
                    pcov_ComFit = pcov_ComFit2
                else:
                    popt_ComFit = popt_ComFit1
                    pcov_ComFit = pcov_ComFit1
                        
            else:
                print('Dataset is nan at: ' + str((x,y)))
                popt_ComFit = np.zeros(np.shape(guesses))
                pcov_ComFit = np.zeros((len(guesses),len(guesses)))
                popt_ComFit[:] = np.nan              
                pcov_ComFit[:] = np.nan
                popt_ComFit1 = popt_ComFit
                popt_ComFit2 = popt_ComFit  
                pcov_ComFit1 = pcov_ComFit
                pcov_ComFit2 = pcov_ComFit 
                Chi2 = np.nan
                red_Chi2 = np.nan
             

            
            # perr_ComFit = np.sqrt(np.diag(pcov_ComFit)) # (6)
            
            
            # (1) calculate the residuals (difference between fit and data 
            #     points). 
            # (2) calculates the normalised residuals, res/err_dep_var, 
            #     in this experiment, the dependant var is the intensity
            # (3) Chi squared statistic, goodness of fit is maximised when 
            #     this is minimised - PHAS0007 UCL
            # (4) K is the number of degrees of freedom, number of variables 
            #     minus number of parameters.
            # (5) A reduced Chi2 statistic close to one is 'good'
            # (6) This calculates the one standard deviation errors of the 
            #     fit parameters since var = sigma^2
        
            #!!! done to here
            # set the paramater names to have the optimised values (or if the
            # optimisation failed, to be NaNs)

            if fit_type == '1Lor':
                [I01, cen1, Gamma1, m, c] = popt_ComFit
                fit_vars = [I01, cen1, Gamma1, m, c] 
            elif fit_type == '2Lor1BWF':
                [I01, cen1, Gamma1, I02, cen2, Gamma2,   
                 I03, cen3, Gamma3, q, m, c] = popt_ComFit
                fit_vars = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                            I03, cen3, Gamma3,  q,  m, c]  
            elif fit_type == '3Lor1BWF':
                [I01, cen1, Gamma1, I02, cen2, Gamma2,   
                 I03, cen3, Gamma3, I04, cen4, Gamma4, 
                 q, m, c] = popt_ComFit
                fit_vars = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                            I03, cen3, Gamma3, I04, cen4, Gamma4, q,  m, c]  
            elif fit_type == '1Lor1BWF':
                [I01, cen1, Gamma1, I02, cen2, Gamma2,
                 q, m, c] = popt_ComFit
                fit_vars = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                            q, m, c]   
            elif fit_type == '1Lor1BWF1Lor':
                [I01, cen1, Gamma1, I02, cen2, Gamma2,   
                 I03, cen3, Gamma3, q, m, c] = popt_ComFit
                fit_vars = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                            I03, cen3, Gamma3, q, m, c]                  
            elif fit_type == '2Lor1BWF1Lor':
                [I01, cen1, Gamma1, I02, cen2, Gamma2,   
                 I03, cen3, Gamma3, I04, cen4, Gamma4, 
                 q, m, c] = popt_ComFit
                fit_vars = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                            I03, cen3, Gamma3, I04, cen4, Gamma4, q, m, c]  
            elif fit_type == '2Lor':
                [I01, cen1, Gamma1, I02, cen2, Gamma2,   
                 m, c] = popt_ComFit
                fit_vars = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                            m, c]  
            elif fit_type == '3Lor':
                [I01, cen1, Gamma1, I02, cen2, Gamma2,
                  I03, cen3, Gamma3, m, c] = popt_ComFit
                fit_vars = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                            I03, cen3, Gamma3,  m, c]  
            elif fit_type == '4Lor':
                [I01, cen1, Gamma1, I02, cen2, Gamma2,
                 I03, cen3, Gamma3, I04, cen4, Gamma4,
                 m, c]  = popt_ComFit
                fit_vars = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                            I03, cen3, Gamma3, I04, cen4, Gamma4,
                            m, c]  
                
            elif fit_type == '5Lor':
                [I01, cen1, Gamma1, I02, cen2, Gamma2,
                  I03, cen3, Gamma3, I04, cen4, Gamma4,
                  I05, cen5, Gamma5, m, c] = popt_ComFit
            
                fit_vars = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                            I03, cen3, Gamma3, I04, cen4, Gamma4,
                            I05, cen5, Gamma5, m, c]  
            
            elif fit_type == '4Lor1BWF':
                [I01, cen1, Gamma1, I02, cen2, Gamma2,
                  I03, cen3, Gamma3, I04, cen4, Gamma4,
                  I05, cen5, Gamma5, q, m, c] = popt_ComFit
            
                fit_vars = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                            I03, cen3, Gamma3, I04, cen4, Gamma4,
                            I05, cen5, Gamma5, q, m, c]  
                
            elif fit_type == '5Lor1BWF':
                [I01, cen1, Gamma1, I02, cen2, Gamma2,
                  I03, cen3, Gamma3, I04, cen4, Gamma4,
                  I05, cen5, Gamma5, I06, cen6, Gamma6, q, m, c] = popt_ComFit
            
                fit_vars = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                            I03, cen3, Gamma3, I04, cen4, Gamma4,
                            I05, cen5, Gamma5, I06, cen6, Gamma6, q, m, c]  
            
            elif fit_type == '6Lor':
                [I01, cen1, Gamma1, I02, cen2, Gamma2,
                  I03, cen3, Gamma3, I04, cen4, Gamma4,
                  I05, cen5, Gamma5, I06, cen6, Gamma6, m, c] = popt_ComFit
            
                fit_vars = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                            I03, cen3, Gamma3, I04, cen4, Gamma4,
                            I05, cen5, Gamma5, I06, cen6, Gamma6, m, c]  
            
            elif fit_type == '7Lor':
                [I01, cen1, Gamma1, I02, cen2, Gamma2,
                  I03, cen3, Gamma3, I04, cen4, Gamma4,
                  I05, cen5, Gamma5, I06, cen6, Gamma6,
                  I07, cen7, Gamma7, m, c] = popt_ComFit
            
                fit_vars = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                            I03, cen3, Gamma3, I04, cen4, Gamma4,
                            I05, cen5, Gamma5, I06, cen6, Gamma6,
                            I07, cen7, Gamma7, m, c]  
                
            elif fit_type == '8Lor':
                [I01, cen1, Gamma1, I02, cen2, Gamma2,
                  I03, cen3, Gamma3, I04, cen4, Gamma4,
                  I05, cen5, Gamma5, I06, cen6, Gamma6,
                  I07, cen7, Gamma7, I08, cen8, Gamma8,
                  m, c] = popt_ComFit
            
                fit_vars = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                            I03, cen3, Gamma3, I04, cen4, Gamma4,
                            I05, cen5, Gamma5, I06, cen6, Gamma6,
                            I07, cen7, Gamma7, I08, cen8, Gamma8,
                            m, c]  
                
            elif fit_type == '9Lor':
                [I01, cen1, Gamma1, I02, cen2, Gamma2,
                  I03, cen3, Gamma3, I04, cen4, Gamma4,
                  I05, cen5, Gamma5, I06, cen6, Gamma6,
                  I07, cen7, Gamma7, I08, cen8, Gamma8,
                  I09, cen9, Gamma9, m, c] = popt_ComFit
                # if Gamma7 < 1:
                #     I07 = 0.001
                fit_vars = [I01, cen1, Gamma1, I02, cen2, Gamma2,
                            I03, cen3, Gamma3, I04, cen4, Gamma4,
                            I05, cen5, Gamma5, I06, cen6, Gamma6,
                            I07, cen7, Gamma7, I08, cen8, Gamma8,
                            I09, cen9, Gamma9, m, c]  
                
            elif fit_type == '2Voigt':
                [amp1, sigma1, I01, cen1, Gamma1,
                 amp2, sigma2, I02, cen2, Gamma2, m, c] = popt_ComFit
                fit_vars = [amp1, sigma1, I01, cen1, Gamma1,
                            amp2, sigma2, I02, cen2, Gamma2, m, c]
            
            elif fit_type == '3Voigt':
                [amp1, sigma1, I01, cen1, Gamma1,
                 amp2, sigma2, I02, cen2, Gamma2,
                 amp3, sigma3, I03, cen3, Gamma3, m, c] = popt_ComFit
                fit_vars = [amp1, sigma1, I01, cen1, Gamma1,
                            amp2, sigma2, I02, cen2, Gamma2,
                            amp3, sigma3, I03, cen3, Gamma3, m, c]
            
            elif fit_type == '4Voigt':
                [amp1, sigma1, I01, cen1, Gamma1,
                 amp2, sigma2, I02, cen2, Gamma2,
                 amp3, sigma3, I03, cen3, Gamma3,
                 amp4, sigma4, I04, cen4, Gamma4, m, c] = popt_ComFit
                fit_vars = [amp1, sigma1, I01, cen1, Gamma1,
                            amp2, sigma2, I02, cen2, Gamma2,
                            amp3, sigma3, I03, cen3, Gamma3, 
                            amp4, sigma4, I04, cen4, Gamma4, m, c]
                
            elif fit_type == '5Voigt':
                [amp1, sigma1, I01, cen1, Gamma1,
                 amp2, sigma2, I02, cen2, Gamma2,
                 amp3, sigma3, I03, cen3, Gamma3,
                 amp4, sigma4, I04, cen4, Gamma4,
                 amp5, sigma5, I05, cen5, Gamma5, m, c] = popt_ComFit
                fit_vars = [amp1, sigma1, I01, cen1, Gamma1,
                            amp2, sigma2, I02, cen2, Gamma2,
                            amp3, sigma3, I03, cen3, Gamma3, 
                            amp4, sigma4, I04, cen4, Gamma4,
                            amp5, sigma5, I05, cen5, Gamma5, m, c]
                
            
            # update dataframe variables, order must match 'p_names' 
            df_var = np.concatenate((fit_vars, np.array([Chi2, red_Chi2])))
           
            # initialise dictionary            
            parameters = dict() 
            # match parameters to dict
            for p in range(len(p_names)):
                parameters[p_names[p]] = df_var[p] 
            
            # populate fit_df with parameters at the map coordinates                       
            fit_df[x,y] = list(parameters.values())
            
    
    # update the overall map library
    fit_data[sample] = fit_df
    
    #----------------------------generate report------------------------------
    

    df = fit_data[sample] # simplify
    x0 = np.zeros(len(p_names)) # initialise
    sigma = np.zeros(len(p_names))
    x0_err = np.zeros(len(p_names))
    counter = 0
    
    chi = df.loc['red_Chi2']
    chi = chi.dropna(axis=0)
    chi = np.array(chi)
    
    for p in p_names: # cycle through parameters

        var = df.loc[p] # for each parameter
        #var = var.drop(labels='Parameters', axis='columns') # drop label
        var = var.dropna(axis=0) # drops NaNs
        var = np.array(var)# convert to array
              
        #---chuck 'bad' data---

        var = var[(np.logical_and(chi>0, chi<thresh_Chi))]
        if len(var) == 0 and counter == 0:
            print('\nWarning: this dataset is outwith thresh_Chi')
            print(p)
            warnings.simplefilter('ignore', category=RuntimeWarning)
        av, stdv = stats.norm.fit(var) 
        var = var[(np.logical_and(var> av - thresh_sig*stdv,
                                  var< av + thresh_sig*stdv))] 
        #----------------------
            
        x0[counter], sigma[counter] = stats.norm.fit(var) # get mean and stdv
        x0_err[counter] = sigma[counter]/np.sqrt(var.shape[0]) # standard error

        counter += 1
        
    report[sample] = pd.DataFrame(index = p_names) # create a dataframe
    report[sample]['x0'] = x0 # save lists to dataframe
    report[sample]['x0_err'] = x0_err
    report[sample]['sigma'] = sigma
    
    #output to excel
    report[sample].to_excel('C:/Users/'+user+'/OneDrive - University College' + \
                        ' London/PhD/Python/Data/Reports/'+ sample + '_' + \
                            fit_type + '_report.xlsx', index = True) 

    if pickle_data == True:
        date = time.strftime("%d%m%y_%H%M")
        name1 = sample + '_raw data_' + date
        name2 = sample + '_clean data_' + date
        name3 = sample + '_fit data_' + date
        names = [name1, name2, name3]
        all_save_data = [raw_data[sample], clean_data[sample], fit_data[sample]]
        paths = list()
        pi = 0
        for name in names:
            path = 'C:/Users/'+user+'/OneDrive - University College London/PhD/' + \
                   'Python/Data/Reports/Pickle Database/' + name + '.pickle'
            data = all_save_data[pi]
            pi+= 1
            with open(path, 'wb') as handle:
                pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    cnt += 1
    
    
    
#-----------------------------------------------------------------------------

# End process timer
end_time = time.perf_counter()
sys.stdout.write('\n')
print("\nScript runtime: %.2f \bs" % (end_time - start_time))
winsound.Beep(frequency, duration)
# last runtime = 360s

#---------------------------------Script End----------------------------------              
