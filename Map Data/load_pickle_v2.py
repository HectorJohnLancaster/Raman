
#=============================================================================
#============================== Pickling Files ===============================
#=============================================================================

# Imports
#import pickle5 as pickle
import pickle
import os

user = os.getlogin()

def load_data(name):
    
    path = 'C:/Users/'+user+'/OneDrive - University College London/PhD/' + \
       'Python/Data/Reports/Pickle Database/' + name + '.pickle'  
       
    with open(path, 'rb') as handle:
        loaded_data = pickle.load(handle)  
        
    return(loaded_data)

# data = ['IE undoped map', 'top1 doped small map1']
# data = ['6 mid burried map1','6 mid doped linemap1']

# data = ['AC kc8 map1']
# data = ['AC og map', 'CB og map', 'Char og map'] 
#data = ['diamond reference map','SC Graphite map','AC og map', 'CB og map', 'Char og map'] 
data = ['SC Graphite map','Graphite S0 MTI map', 'AC og map', 'CB og map', 'Char og map']
#data = ['AC og map', 'AC kc36 map', 'AC kc30 map', 'AC kc24 map', 'AC kc10 map', 'AC kc8 map', 'AC kc4 map']
# data = ['AC og map', 'AC kc36 map', 'AC kc36 discharged map', 
#         'AC kc30 map', 'AC kc30 discharged map',
#         'AC kc24 map', 'AC kc24 discharged map',
#         'AC kc10 map', 'AC kc10 discharged map',
#         'AC kc8 map', 'AC kc8 discharged map',
#         'AC kc4 map', 'AC kc4 discharged map']
# data = ['AC og map', 'AC kc10 map', 'AC kc10 discharged map']
# data = ['AC og map', 'AC kc36 map', 'AC kc30 map', 'AC kc24 map', 'AC kc10 map', 'AC kc8 map', 'AC kc4 map',
#         'CB og map', 'CB kc24 redo map', 'CB kc10 map', 'CB kc8 map',
#         'Char og map', 'Char kc24 map', 'Char kc10 map', 'Char kc8 map'] 
# data = ['Char og map', 'Char kc24 map', 'Char kc10 map', 'Char kc8 map']
# data = ['Char og map', 'Char kc24 map', 'Char kc24 discharged map', 'Char kc10 map',
#         'Char kc10 discharged map', 'Char kc8 map', 'Char kc8 discharged map']
# data = ['CB og map', 'CB kc24 redo map', 'CB kc24 redo discharged map', 'CB kc10 map', 
#         'CB kc10 discharged map', 'CB kc8 map', 'CB kc8 discharged map']
#data = ['AC og map', 'AC kc10 map', 'AC kc10 discharged map', 'AC NaC10 map', 'AC NaC10 discharged map']

raw_data = dict()
clean_data = dict()
fit_data = dict()
sample_names = list()
for i in data:
    raw_data[i] = load_data(i+'_nas_raw data') # nas = no added signal
    clean_data[i] = load_data(i+'_nas_clean data')
    fit_data[i] = load_data(i+'_nas_fit data')
    # raw_data[i] = load_data(i+'_raw data')
    # clean_data[i] = load_data(i+'_clean data')
    # fit_data[i] = load_data(i+'_fit data')
    sample_names.append(i)

#=============================================================================
#=============================================================================