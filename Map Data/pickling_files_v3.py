
#=============================================================================
#============================== Pickling Files ===============================
#=============================================================================

#----------------------------------Preamble-----------------------------------

# File description: 
# Saving a dictionary as a pickle

# Imports
import pickle
        
#---------------------------------User Input----------------------------------

# dictionary to save
#sample = sample_names[1]
dsave = ['raw data', 'clean data', 'fit data']


#----------------------------------The Meat-----------------------------------
for sample in sample_names:
    def yes_or_no(question):
        while "the answer is invalid":
            reply = str(input(question+' (y/n): ')).lower().strip()
            if reply[0] == 'y':
                answer = True
                return answer
            if reply[0] == 'n':
                answer = False
                return answer
    
    for data in dsave:   
        #------     
        name = sample + '_nas_' + data
        path = 'C:/Users/'+user+'/OneDrive - University College London/PhD/' + \
               'Python/Data/Reports/Pickle Database/' + name + '.pickle'
        if data == 'raw data':       
            save_dict = raw_data[sample]
        elif data == 'clean data':
            save_dict = clean_data[sample]
        elif data == 'fit data':
            save_dict = fit_data[sample]
        #------     
        try:
            with open(path, 'xb') as handle:
                pickle.dump(save_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
            print("\nFile saved")
            
        except FileExistsError:
            answer = yes_or_no('This file already exists, are you sure you want to overwrite?')
            if answer == True:
                with open(path, 'rb') as handle:
                    loaded_data = pickle.load(handle)
                print("\nFile saved")
            if answer == False:
                print("\nSave aborted")
            

#=============================================================================
#=============================================================================