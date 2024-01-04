This series of code is capable of processing Raman point data from a txt file in the format given in "example_point_data.txt". 

**The code must be run in the following order:**
1. **pointdata_process**  
   info: imports point data from txt file and stores in accessible format and allows for normalisation, rebinning, cosmic-ray removal and fluorescent background removal  
   requirements: if rebinning requires "xwindow" in directory; if fluorescent background removal requires "rolling_ball_single" in directory
2. **pointdata_fit**  
   info: fits point data with defined peaks  
   requirements: requires "initial_guesses" in directory  

Optional code:  
1. rolling_ball_single_checker (allows for a given rolling ball background removed map spectra to be investigated and visualised)
