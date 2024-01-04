
#-----------------------------import modules---------------------------------

import time
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import initial_guesses # Must be in the working directory

def spectrum_interval(x,y, lb, ub):
    lowerx = np.where(x>lb)
    upperx = np.where(x<ub)
    interx = np.intersect1d(lowerx,upperx)
    x = x[interx]  
    y = y[interx]
    return(x,y)


#====== From Map Scan =======
coord = (80,0)#(-0.38, 1.6)#(-2,8)#(-36,-8)#(50,55)#
sample = sample_names[0]
xs = clean_data[sample][coord][2,:]#xspec
ys = clean_data[sample][coord][3,:]# + 10000
#============================

# #===== From Point Scan ======
# coord = sample_names[1]
# xs = clean_data[coord][:,0]
# ys = clean_data[coord][:,1] #+ 100000
# #============================

#======== Other Crap ========
# #coord = 'graphite'
# # coord = (-12,8)
xs, ys = spectrum_interval(xs, ys, 1150, 3200)

#plt.ylim(0,310)
# xs = sliced_data[sample][coord][:,2]#xspec
# ys = sliced_data[sample][coord][:,3] + 10000

y2 = np.mean(ys[:10])
y1 = np.mean(ys[-10:])
x2 = xs[0]
x1 = xs[-1]
m = (y2-y1)/(x2-x1)     
c = float(y2 - m*x2)   
xtemp = np.arange(1150,3200,10)

# plt.plot(xs,ys)
# plt.plot(xtemp, xtemp*m+c)

# xs = xs2
# ys = t2 + 10

# coord =''
# xs = xsafe
# ys = ysafe #+ 10000

# xs = x
# ys = y
#============================

#xs, ys = spectrum_interval(xs, ys, 950, 1850)


#--------------------------------user inputs----------------------------------

# Choose which data to chuck, below this threshold of reduced Chi2 (thresh_Chi)
# and fit_var standard deviations (thresh_sig) is kept.
thresh_Chi = 20
thresh_sig = 6

fit_type_list = ['2Lor', '3Lor', '4Lor', '5Lor', '6Lor', '7Lor', '8Lor', '9Lor',
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

fit_type = '1Lor1BWF1Lor'#'2Lor1BWF'#'3Lor1BWF'#
material_type = 'SC Graphite'#'Active Carbon'#'graphite'#'AC kc8'#'AC kc8'#'MnOx4'#material_list[8]


#----------------------------start process timer------------------------------

start_time = time.perf_counter()


#-----------------------------define functions--------------------------------


#_____________________________individual curves_______________________________

# See cauchy distribution function at:
# https://www.itl.nist.gov/div898/handbook/eda/section3/eda3663.htm 
# Here, the  modified version with the max intensity (I0) replaces the 
# the parameter (1/s*pi).
def Lorentzian(x, I0, cen, Gamma, m, c):
    return I0*(Gamma/2)**2/((x-cen)**2+(Gamma/2)**2) + (m * x + c)

# Function as described in SI of DOI: 10.1103/PhysRevB.84.241404    
def BWF(x, I0, cen, Gamma, q, m, c):
    numerator = (1+(x-cen)/(q*(Gamma/2)))**2
    denominator = 1+((x-cen)/(Gamma/2))**2
    return I0*(numerator/denominator) + (m * x + c)

def Voigt(x, amp, sigma, I0, cen, Gamma, m, c):
    return (amp*(1/(sigma*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen)**2)/(2*sigma**2)))) +\
            ((I0*(Gamma/2)**2/((x-cen)**2+(Gamma/2)**2)) ) + (m * x + c)

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


#xs, ys = spectrum_interval(xs, ys, 1000, 1800)
#ys = ys/max(ys)
#xs, ys = xwindow(xs, ys, 1)

#_____________________________combined function_______________________________

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

     # Choose parameter names 'p_names' for the dataframe, order must match  
     # the order that the parameters are storred upon entry.
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
    
    
#--------------------------get & print fitting data---------------------------


noise = 0
            
if material_type == 'Alex Diamond':
    xs = xs[ys != 0] # delete any zeros
    ys = ys[ys != 0]
                
            
ys_err = np.sqrt(abs(ys)) # (1)
ys_err += 1*10**-20 # (2)
            
# (1) Uncertainty in intensity values, since photons hitting a 
#     detector follow a Poisson distribution.
# (2) Adds a tiny uncertainty to avoid dividing by zero later on


#---------------------find BWF/Lorentzian fit---------------------   


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
    fig, ax = plt.subplots()
    ax.plot(xs, ys, label= 'material: ' + \
            sample)
    ax.legend()

    if second_pass == True:
        checker = True
        print('RuntimeError, optimisation only valid on first pass')
        popt_ComFit2 = np.zeros(np.shape(guesses))
        pcov_ComFit2 = np.zeros((len(guesses),len(guesses)))
        popt_ComFit2[:] = np.nan              
        pcov_ComFit2[:] = np.nan

    else:
        print('optimisation failed')
        popt_ComFit = np.zeros(np.shape(guesses))
        pcov_ComFit = np.zeros((len(guesses),len(guesses)))
        popt_ComFit[:] = np.nan              
        pcov_ComFit[:] = np.nan
        popt_ComFit1 = popt_ComFit
        popt_ComFit2 = popt_ComFit  
        pcov_ComFit1 = pcov_ComFit
        pcov_ComFit2 = pcov_ComFit 

except ValueError:
        print('ValueError, optimisation failed')
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
        print('Data too noisy, set to NaN')
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

if sred_Chi2 < red_Chi2:#15:#+0.05:
    popt_ComFit = popt_ComFit2
    pcov_ComFit = pcov_ComFit2
else:
    popt_ComFit = popt_ComFit1
    pcov_ComFit = pcov_ComFit1

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

#--------------------------store paramaters-----------------------


# ratio of intensity maxima, I0D/I0G
p1 = [I01, cen1, Gamma1]
p2 = [I02, cen2, Gamma2]
Iratio = p1[0]/p2[0]

# integrate individual curves to get their area, here the bounds
# are taken from the range of the scan


    
# update dataframe variables, order must match 'p_names' 
df_var = np.concatenate((fit_vars, np.array([Chi2, 
                                              red_Chi2])))
   


# initialise dictionary            
parameters = dict() 
# match parameters to dict
for p in range(len(p_names)):
    parameters[p_names[p]] = df_var[p] 

#q=-50
fig, ax = plt.subplots(figsize=(3.6,2.7))
        
background = m * xs + c   
background = 0 * xs + 0
ys = ys - background 

nBWF = list()
nLor = list()
iBWF = list()
iLor = list()
icheck = 0
for i in range(fit_type.count('BWF')):
    iBWF.append(fit_type.index('BWF', icheck) - 1)
    icheck += 3
    nBWF.append(int(fit_type[iBWF[i]]))
icheck = 0
for i in range(fit_type.count('Lor')):
    iLor.append(fit_type.index('Lor', icheck) - 1)
    icheck += 3
    nLor.append(int(fit_type[iLor[i]]))

fv = fit_vars
# draft code for plotting any alternating combo of LorBWFLor or BWFLorBWF
if len(iBWF) == 0: # if no BWF
    for n in range(nLor[0]): # for each Lorentzian
        ln = np.array([fv[0+3*n], fv[1+3*n], fv[2+3*n], m, c]) # get fit vars
        lorn = Lorentzian(xs, *ln) - background # get y data
        ax.plot(xs, lorn, color="grey", alpha=0.5, linewidth=0.5) # plot
if len(iLor) == 0: # if no Lorentzians
    for n in range(nBWF[0]):
        bwfn = np.array([fv[0+3*n], fv[1+3*n], fv[2+3*n], q, m, c])
        BWFn = BWF(xs, *bwfn) - background
        ax.plot(xs, BWFn, color="grey", alpha=0.5, linewidth=0.5)
try: # if there is a combo of both Lorentzians and BWF lineshapes
    if iLor[0] < iBWF[0]: # if lor first
        for n in range(nLor[0]): 
            ln = np.array([fv[0+3*n], fv[1+3*n], fv[2+3*n], m, c]) 
            lorn = Lorentzian(xs, *ln) - background 
            ax.plot(xs, lorn, color="grey", alpha=0.5, linewidth=0.5) 
            Lor_first = True
    else: # else if bwf first
        for n in range(nBWF[0]):
            bwfn = np.array([fv[0+3*n], fv[1+3*n], fv[2+3*n], q, m, c])
            BWFn = BWF(xs, *bwfn) - background
            ax.plot(xs, BWFn, color="grey", alpha=0.5, linewidth=0.5)
            Lor_first = False
    if Lor_first == True: # if lor first, do bwf second
        nprime = nLor[0]
        for n in range(nBWF[0]):
            bwfn = np.array([fv[0+3*(n+nprime)], fv[1+3*(n+nprime)], fv[2+3*(n+nprime)], q, m, c])
            BWFn = BWF(xs, *bwfn) - background
            ax.plot(xs, BWFn, color="grey", alpha=0.5, linewidth=0.5)
            BWF_second = True   
    else: # otherwise if bwf first, do lor second
        nprime = nBWF[0]
        for n in range(nLor[0]):
            ln = np.array([fv[0+3*(n+nprime)], fv[1+3*(n+nprime)], fv[2+3*(n+nprime)], m, c])
            lorn = Lorentzian(xs, *ln) - background
            ax.plot(xs, lorn, color="grey", alpha=0.5, linewidth=0.5)
            BWF_second = False
    if BWF_second == True: # if bwf second, do lor next
        nprime = nLor[0]+nBWF[0]
        for n in range(nLor[1]):
            ln = np.array([fv[0+3*(n+nprime)], fv[1+3*(n+nprime)], fv[2+3*(n+nprime)], m, c])
            lorn = Lorentzian(xs, *ln) - background
            ax.plot(xs, lorn, color="grey", alpha=0.5, linewidth=0.5)
    else: # else if lor second, do bwf next
        nprime = nLor[0]+nBWF[0]
        for n in range(nBWF[1]):
            bwfn = np.array([fv[0+3*(n+nprime)], fv[1+3*(n+nprime)], fv[2+3*(n+nprime)], q, m, c])
            BWFn = BWF(xs, *bwfn) - background
            ax.plot(xs, BWFn, color="grey", alpha=0.5, linewidth=0.5) 
except IndexError: # if index error due to exceeding number of BWF/Lors, pass
    pass
    
comfit = ComFit(xs, *fit_vars) - background

#---

#xs, ys = spectrum_interval(xs, ys, 800, 1800)
ax.plot(xs, ys, linestyle='-', marker = 'o', markersize = 0.8, color="k")#linestyle='', marker = 'o', mew=0, markersize =1.2, mfc="m", label= 'Raw Data')
ax.plot(xs, comfit, linestyle='-', color="r", label= 'fit', alpha=0.6, linewidth=1)


# ax.text(0.02,0.7, 'cen1: %.f \ncen2: %.f \ncen3: %.f \nred$\chi$$^{2}$: %.3f' % (cen1, cen2, cen3, red_Chi2),
#         transform=ax.transAxes, color='r')
# ax.text(0.05,0.7, 'cenD: %.f \ncenG: %.f \nFano: %.2f \n$\chi$$^{2}$: %.2f' % (cen1, cen3, -1/q, Chi2),
#         transform=ax.transAxes, color='r')

#---
#ax.legend()
ax.set_xlabel('Raman shift (cm$^{-1}$)')
ax.set_ylabel('Intensity (arb. units)')
ax.set_yticks([])
#plt.title('Example Fit: ' + str(coord))
fig.tight_layout()

ax.set_xlim(2500,2800)
#ax.set_ylim(-40,200)

#print(I03/I01)

fig.savefig('C:/Users/'+user+'/OneDrive - University College London' + \
            '/PhD/Python/Data/Figures/Gyen eg fit'+str(coord)+'.svg', dpi=300,
            transparent = True)
#-----------------------------------------------------------------------------


#---------------------------------Script End----------------------------------              
