#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import necesarry bibs
import numpy as np
import scipy.io as scio
import matplotlib.pyplot as plt
from PIL import Image
from IPython import get_ipython


def pupil_resize(desired_size, p_custom_pupil):
    """
    Created on Tue Dec 18 15:56:38 2018

    @author: hadorn
    
    This fucntion rescales a pupil, based on the given .mat File
    
    If the resulting plot is still shown in consol, change 
    Tools >> Preferences >> IPython console >> Graphics --> Automatic
    
    Inputs: 
        desired_size : (double) : Scaling of the resulting pupil 
        p_custom_pupil : (custom pupil obejct) : Parameter object of custom pupil  
    
    Outputs:
        resized_pupil: Scaled pupil as array
    """ 
    
    if p_custom_pupil.custom_pupil_path[-3:] != 'mat':
        
        raise TypeError("custom pupil data must be a .mat file.")
    
    else:
    
        pupil_dic = scio.loadmat(p_custom_pupil.custom_pupil_path) # load pupil data from .mat file 
    
        pupil = pupil_dic[p_custom_pupil.custom_pupil_name] # get pupil from pupil_dic
    
        scale_factor = desired_size/len(pupil) # determine scale factor
            
        scaled_n = int(len(pupil)*scale_factor) # get length of pupil and determine new scaled length
                
        resized_pupil = np.array(Image.fromarray(pupil).resize((scaled_n,scaled_n))) # rescale pupil
               
        if p_custom_pupil.show_pupil == True: # check if figure is needed
                
            # plot pupil in external window, if it does not, reading the header might be useful
            get_ipython().run_line_magic('matplotlib','qt')
                                
            plt.imshow(resized_pupil)
            
            
            return resized_pupil # return resized pupil
            
        else:
                
            return resized_pupil # return resized pupil