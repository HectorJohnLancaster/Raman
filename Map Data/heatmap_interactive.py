# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 10:17:22 2022

@author: hlancaster
"""

import holoviews as hv
import hvplot.pandas # despite error, this is used 
hv.extension('bokeh')  # to generate interactive plots 

sample = sample_names[0]

hm = mat_grid[sample].hvplot.heatmap(cmap='inferno', flip_yaxis=True,  
                         xaxis=None, yaxis=None,
                         frame_width=300, frame_height=300,
                         clabel = 'cbar')
hv.save(hm, 'C:/Users/'+user+'/OneDrive - University College London' + \
            '/PhD/Python/Data/Figures/heatmap.html')