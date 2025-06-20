This series of code is capable of processing Raman mapping data from a txt file in the format given in "example_map_data.txt". 

**The code must be run in the following order:**
1. **map_process**  
   info: imports map data from txt file and stores in accessible format and allows for normalisation, rebinning, cosmic-ray removal and fluorescent background removal  
   requirements: if rebinning requires "xwindow" in directory; if fluorescent background removal requires "rolling_ball_map" in directory
3. **map_fit**  
   info: fits map data with defined peaks  
   requirements: requires "initial_guesses" in directory
5. **heatmap**  
   info: produces a heatmap of specified fit parameters form the fitted mapping data  
6. **histogram**  
   info: produces a histogram based on the data included in the heatmap  

Optional code: 
1. area (calculates integrated area ratios and peak intensity ratios from fitted data, to be run after "map_fit")
2. heatmap_interactive (produces an interactive heatmap, loaded in a browser window, allowing for pixels to be examined for their x,y coordinates)
3. pickling_files (allows processed and fitted data from "map_process" and "map_fit" to be stored for future use
4. load_pickle (allows for pickled data from "pickling_files" to be loaded in for use)
5. rolling_ball_single_checker (allows for a given rolling ball background removed map spectra to be investigated and visualised)
