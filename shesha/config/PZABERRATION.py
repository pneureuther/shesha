"""
Param_z_aberration class definition

Parameters for Zernike Aberration Utility

@author: neureuther
"""

import numpy as np


#################################################
# P-Class (parametres) Zernike Aberration
#################################################
class Param_z_aberration:
    
    #%%
    #################################################
    # DEFINITION
    #################################################
    def __init__(self):
        """
        Parameters to configure the additional Zernike aberrations
        include_zab     if true: Zernike aberrations are added to pupil phase 
                        in COMPASS simulations (bool)
        file_dir        directory of mat-file containing the time series of 
                        Zernike coefficients (string)
        file_name       file-name of mat-file containing the time series of 
                        Zernike coefficients (string)
        num_zpol        Number of Zernike polynomials (Noll convention) 
                        included in simulation (integer)
        zcube_spup      Cube of Zernike modes for spupil (scientific path/target)
                        (np.ndarray[ndim=3, dtype=np.float64)
        zcube_mpup      Cube of Zernike modes for mpupil (analytic path/WFS)
                        (np.ndarray[ndim=3, dtype=np.float64)    
        coeff_wfs       Timeseries of Zernike coefficients without time stamps 
                        for wfs path (np.ndarray[ndim=2, dtype=np.float64)
        coeff_sci       Timeseries of Zernike coefficients without time stamps 
                        for science path (np.ndarray[ndim=2, dtype=np.float64)
        time_series     Time stamps (in seconds) of timeseries of Zernike
                        coefficients (np.ndarray[ndim=1, dtype=np.float64)
        include_path    Variable to define in which paths the aberrations are
                        taken into account (integer)
                        0 = not included in science and analytick path
                        1 = included only in science path (target)
                        2 = included only in analytic path (WFS)
                        3 = included in science and analytic path
        mat_vers        Version of the mat-file containing the timeseries of 
                        Zernike coefficients (string)
        step            Time steps (in seconds) of time_series (float64)
        var_name_c_wfs  Variable-name of coeff_wfs in mat-file (string)
        var_name_c_sci  Variable-name of coeff_sci in mat-file (string)
        var_name_time   Variable-name of time_series in mat-file (string)
        diam_data       Diameter (in meter) of the aberrations stored in the 
                        mat-file (float)
        pup_diam        Diameter (in meter) of the telescope pupil to which the
                        the aberrations are truncated (float)
                        -2. = ignore different diameter and add complete Zernike
                              modes to the simulation
                        -1. = inherit diameter from p_tel.diam
                        >0. = user specified telescope pupil diameter
        dec             decimation of iteration time to aberration time steps
                        (= sim.config.p_loop.it_time / step) (int)
        phase_index     Index of phasescreen aberrations (int)
        """
        self.__include_zab = False
        self.__file_dir = "default"
        self.__file_name = "input.mat"
        self.__num_zpol = 0
        self.__zcube_spup = None
        self.__zcube_mpup = None
        self.__coeff_wfs = None
        self.__coeff_sci = None
        self.__time_series = None
        self.__include_path = 0
        self.__mat_vers = "v7.3"
        self.__step = 0.0
        self.__var_name_c_wfs = "coeff_wfs"
        self.__var_name_c_sci = "coeff_science"
        self.__var_name_time = "time"
        self.__diam_data = 0.0
        self.__pup_diam = -2.0
        self.__dec = 0
        self.__phase_index = 0
    
    #%%
    #################################################
    # SET-COMMANDS
    #################################################
    def set_include_zab(self, include_bool):
        """ Set if Zernike aberration will be included
        
        :parameters:
            include_bool: (bool) : bool determining if Zernike aberration will
                                   be included
        """
        if type(include_bool) != bool:
            raise TypeError("include_bool must be a boolean (True or False)")
        
        self.__include_zab = include_bool
        
    include_zab = property(lambda x: x.__include_zab, set_include_zab)
    
    #%%
    def set_file_dir(self, f_dir):
        """ Set the directory of mat-file containing the time series of Zernike
            coefficients
        
        :parameters:
            f_dir: (string) : directory of mat-file containing the time series
                              of Zernike coefficients
        """
        if type(f_dir) != str:
            raise TypeError("file_dir must be a string")
        
        self.__file_dir = f_dir
    
    file_dir = property(lambda x: x.__file_dir, set_file_dir)
    
    #%%
    def set_file_name(self, f_name):
        """ Set the file-name of mat-file containing the time series of Zernike
            coefficients
        
        :parameters:
            f_name: (string) : file-name of mat-file containing the time series
                               of Zernike coefficients
        """
        if type(f_name) != str:
            raise TypeError("file_name must be a string")
        
        self.__file_name = f_name
    
    file_name = property(lambda x: x.__file_name, set_file_name)
    
    #%%
    def set_num_zpol(self, num):
        """ Set the number of Zernike polynomials included in simulation
        
        :parameters:
            num: (int) : number of Zernike polynomials included in simulation
        """
        if type(num) != int:
            raise TypeError("num_zpol must be a POSITIVE integer")
        if num <= 0:
            raise TypeError("num_zpol must be a POSITIVE integer")
        
        self.__num_zpol = num
    
    num_zpol = property(lambda x: x.__num_zpol, set_num_zpol)
    
    #%%
    def set_zcube_spup(self, cube):
        """ Set the cube of Zernike modes for spupil (scientific path/target)
        
        :parameters:
            cube: (np.ndarray[ndim=3, dtype=np.float64) : Cube of Zernike modes
                                      for spupil
        """
        if type(cube) != np.ndarray:
            raise TypeError("zcube_spup must be a numpy.ndarray")
        if len(cube.shape) != 3:
            raise TypeError("zcube_spup must be a ndarray with exactly 3 dimensions")
        if len(cube[:,0,0]) != self.__num_zpol:
            raise TypeError("zcube_spup must be a ndarray whose first dimension has the same number of entries as num_zpol")
        if type(cube[0,0,0]) != np.float64:
            raise TypeError("zcube_spup must be a ndarray containig np.float64 numbers")
        
        self.__zcube_spup = cube
    
    zcube_spup = property(lambda x: x.__zcube_spup, set_zcube_spup)
    
    #%%
    def set_zcube_mpup(self, cube):
        """ Set the cube of Zernike modes for mpupil (analytic path/wfs)
        
        :parameters:
            cube: (np.ndarray[ndim=3, dtype=np.float64) : Cube of Zernike modes
                                      for mpupil
        """
        if type(cube) != np.ndarray:
            raise TypeError("zcube_mpup must be a numpy.ndarray")
        if len(cube.shape) != 3:
            raise TypeError("zcube_mpup must be a ndarray with exactly 3 dimensions")
        if len(cube[:,0,0]) != self.__num_zpol:
            raise TypeError("zcube_mpup must be a ndarray whose first dimension has the same number of entries as num_zpol")
        if type(cube[0,0,0]) != np.float64:
            raise TypeError("zcube_mpup must be a ndarray containig np.float64 numbers")
        
        self.__zcube_mpup = cube
    
    zcube_mpup = property(lambda x: x.__zcube_mpup, set_zcube_mpup)
    
    #%%
    def set_coeff_wfs(self, coeff):
        """ Set the timeseries of Zernike coefficients without time stamps for 
        the wfs path
        
        :parameters:
            coeff: (np.ndarray[ndim=2, dtype=np.float64) : timeseries of Zernike
                                       coefficients without time
        """
        if type(coeff) != np.ndarray:
            raise TypeError("coeff_wfs must be a numpy.ndarray")
        if len(coeff.shape) != 2:
            raise TypeError("coeff_wfs must be a ndarray with exactly 2 dimensions")
        if len(coeff[0,:]) != self.__num_zpol:
            raise TypeError("coeff_wfs must be a ndarray whose second dimension has the same number of entries as num_zpol")
        if type(coeff[0,0]) != np.float64:
            raise TypeError("coeff_wfs must be a ndarray containig np.float64 numbers")
        
        self.__coeff_wfs = coeff
    
    coeff_wfs = property(lambda x: x.__coeff_wfs, set_coeff_wfs)
    
    #%%
    def set_coeff_sci(self, coeff):
        """ Set the timeseries of Zernike coefficients without time stamps for 
        the science path
        
        :parameters:
            coeff: (np.ndarray[ndim=2, dtype=np.float64) : timeseries of Zernike
                                       coefficients without time
        """
        if type(coeff) != np.ndarray:
            raise TypeError("coeff_sci must be a numpy.ndarray")
        if len(coeff.shape) != 2:
            raise TypeError("coeff_sci must be a ndarray with exactly 2 dimensions")
        if len(coeff[0,:]) != self.__num_zpol:
            raise TypeError("coeff_sci must be a ndarray whose second dimension has the same number of entries as num_zpol")
        if len(coeff[:,0]) != len(self.coeff_wfs[:,0]):
            raise TypeError("coeff_sci must be a ndarray whose first dimension has the same length as the first dimension of coeff_wfs")
        if type(coeff[0,0]) != np.float64:
            raise TypeError("coeff_sci must be a ndarray containig np.float64 numbers")
        
        self.__coeff_sci = coeff
    
    coeff_sci = property(lambda x: x.__coeff_sci, set_coeff_sci)
    
    #%%
    def set_time_series(self, time):
        """ Set the time stamps (in seconds) of timeseries of Zernike coefficients
        
        :parameters:
            time: (np.ndarray[ndim=1, dtype=np.float64) : timeseries of Zernike
                                      coefficients without time
        """
        if type(time) != np.ndarray:
            raise TypeError("time_series must be a numpy.ndarray")
        if len(time.shape) != 1:
            raise TypeError("time_series must be a ndarray with exactly 1 dimensions")
        if time.size != len(self.coeff_wfs[:,0]):
            raise TypeError("time_series must be a ndarray whose first dimension has the same length as the first dimension of coeff_wfs")
        if type(time[0]) != np.float64:
            raise TypeError("time_series must be a ndarray containig np.float64 numbers")
        
        self.__time_series = time
    
    time_series = property(lambda x: x.__time_series, set_time_series)

    #%%
    def set_include_path(self, path):
        """ Set the variable controlling in which paths the aberrations are 
            taken into account
        
        :parameters:
            path: (int) : variable controlling in which paths the aberrations
                          are taken into account
        """
        if type(path) != int:
            raise TypeError("include_path must be an integer")
        if not path in [0, 1, 2, 3]:
            raise TypeError("include_path must be 0 or 1 or 2 or 3")
        
        self.__include_path = path
    
    include_path = property(lambda x: x.__include_path, set_include_path)
    
    #%%
    def set_mat_vers(self, ver):
        """ Set the version of the mat-file containing the timeseries of 
            Zernike coefficients
        
        :parameters:
            ver: (string) : version of the mat-file containing the timeseries 
                            of Zernike coefficients
        """
        if type(ver) != str:
            raise TypeError("mat_vers must be a string")
        if not ver in ["v4", "v6", "v7", "v7.3"]:
            raise TypeError("mat_vers must be v4 or v6 or 7 or v7.3")
        
        self.__mat_vers = ver
        
    mat_vers = property(lambda x: x.__mat_vers, set_mat_vers)
    
    #%%
    def set_step(self, t):
        """ Set the time step (in seconds) of time_series
        
        :parameters:
            t: (np.float64) : time step (in seconds) of time_series
        """
        if type(t) != np.float64:
            raise TypeError("step must be a POSITIVE numpy.float64")
        if t <= 0:
            raise TypeError("step must be a POSITIVE numpy.float64")
        
        self.__step = t
    
    step = property(lambda x: x.__step, set_step)
    
    #%% 
    def set_var_name_c_wfs(self, name):
        """ Set the variable-name of coeff_wfs in mat-file
        
        :parameters:
            name: (string) : variable-name of coeff_wfs in mat-file
        """
        if type(name) != str:
            raise TypeError("var_name_c_wfs must be a string")
        
        self.__var_name_c_wfs = name
    
    var_name_c_wfs = property(lambda x: x.__var_name_c_wfs, set_var_name_c_wfs)
    
    #%% 
    def set_var_name_c_sci(self, name):
        """ Set the variable-name of coeff_sci in mat-file
        
        :parameters:
            name: (string) : variable-name of coeff_sci in mat-file
        """
        if type(name) != str:
            raise TypeError("var_name_c_sci must be a string")
        
        self.__var_name_c_sci = name
    
    var_name_c_sci = property(lambda x: x.__var_name_c_sci, set_var_name_c_sci)
    
    #%%
    def set_var_name_time(self, name):
        """ Set the variable-name of time_series in mat-file
        
        :parameters:
            name: (string) : variable-name of time_series in mat-file
        """
        if type(name) != str:
            raise TypeError("var_name_time must be a string")
        
        self.__var_name_time = name
    
    var_name_time = property(lambda x: x.__var_name_time, set_var_name_time)
    
    #%%
    def set_diam_data(self, diam):
        """ Set the diameter of aberrations (in meter) stored in mat-file
        
        :parameters:
            diam: (float) : diameter of aberrations (in meter)
        """
        if type(diam) != float:
            raise TypeError("diam_data must be a POSITIVE float")
        if diam <= 0:
            raise TypeError("diam_data must be a POSITIVE float")
        
        self.__diam_data = np.float64(diam)
        
    diam_data = property(lambda x: x.__diam_data, set_diam_data)
    
    #%%
    def set_pup_diam(self, diam):
        """ Set the diameter (in meter) of the telescope pupil
        
        :parameters:
            diam: (float) : diameter (in meter) of the telescope pupil
        """
        if (type(diam) != float) and (type(diam) != np.float64):
            raise TypeError("diam_data must be a float")
        if (diam <= 0):
            if (diam != -1.0) and (diam != -2.0):
                raise TypeError("pup_diam must be a float >0. or =-1. or =-2.")
        
        self.__pup_diam = np.float64(diam)
        
    pup_diam = property(lambda x: x.__pup_diam, set_pup_diam)
    
    #%%
    def set_dec(self, D):
        """ Set the decimation of iteration time to aberration time steps
        
        :parameters:
            D: (integer) : decimation of iteration time to aberration time steps
        """
        
        if (type(D) != int) and (type(D) != np.int64):
            raise TypeError("dec must be a POSITIVE integer.")
        if (D <= 0):
            raise TypeError("dec must be a POSITIVE integer.")
        
        self.__dec = int(D)
    
    dec = property(lambda x: x.__dec, set_dec)
    
    #%%
    def set_phase_index(self, idx):
        """ Set index of phasescreen aberrations
        
        :parameters:
            idx: (integer) : Index of phasescreen aberrations
        """
        if (type(idx) != int) and (type(idx) != np.int64):
            raise TypeError("phase_index must be a NON-NEGATIVE integer.")
        if idx < 0:
            raise TypeError("phase_index must be a NON-NEGATIVE integer.")
        
        self.__phase_index = int(idx)
    
    phase_index = property(lambda x: x.__phase_index, set_phase_index)