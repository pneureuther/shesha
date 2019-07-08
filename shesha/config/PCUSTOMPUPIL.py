#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 14 14:08:32 2019

@author: hadorn

This parameter class enables to get a customized pupil from another directory

"""

class Param_custom_pupil:

    def __init__(self):
        """ import pupil"""
        self.__import_pupil = False
        
        """ pupil name in .mat file"""
        self.__custom_pupil_name = None
        
        """ folder of custom pupil file"""
        self.__custom_pupil_path = None
        
        """ boolean variable for showing the pupil"""
        self.__show_pupil = False
        
    def set_import_pupil(self, imp_pupil, p_tel):
        """ determine if a image of the pupil is wanted

        :param img: (boolean) : show pupil
        :param p_tel: (telescope object) : object of telescope class"""
        
        if type(imp_pupil) != bool:
            raise TypeError("import_pupil must be a boolean.")

        self.__import_pupil = imp_pupil
        
        # set custom pupil aperture type
        if imp_pupil == True:
            
            p_tel.set_type_ap(b'Custom-Pupil')  
    
    import_pupil = property(lambda x: x.__import_pupil, set_import_pupil)    
    
    def set_custom_pupil_name(self, name):
        """ Set the pupil name in .mat file

        :param name: (string) : name of pupil in .mat file """
        
        if type(name) != str:
            raise TypeError("custom_pupil_name must be a string.")
        
        self.__custom_pupil_name = name

    custom_pupil_name = property(lambda x: x.__custom_pupil_name, set_custom_pupil_name)

    def set_custom_pupil_path(self, path):
        """ Set path of the .mat file

        :param path: (string) : path to directory of .mat file """
        
        if type(path) != str:
            raise TypeError("custom_pupil_path must be a string.")
        
        self.__custom_pupil_path = path
    
    custom_pupil_path = property(lambda x: x.__custom_pupil_path, set_custom_pupil_path)   
    
    def set_show_pupil(self, img):
        """ determine if a image of the pupil is wanted

        :param img: (boolean) : show pupil"""
        
        if type(img) != bool:
            raise TypeError("img must be a boolean.")
        
        self.__show_pupil = img
    
    show_pupil = property(lambda x: x.__show_pupil, set_show_pupil)