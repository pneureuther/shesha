"""
Zernike Aberration Utility

Introduce custom Zernike aberrations to simulation

@author: neureuther
"""

import os
import numpy as np
import scipy.io as scio
import hdf5storage
import math

#%%
def calc_zcube(n_zpol, grid_p):
    
    """
    Calculates a Cube of Zernike polynomials. The polynomials can be accessed
    with Z[k,:,:] and are sorted according to Noll convention.
    Z[0,:,:] = pistion
    Z[1,:,:] = tip
    Z[2,:,:] = tilt
    [...]
    
    :parameters:
        n_zpol: (integer):  Number of Zernike polynomials
        grid_p: (integer):   Number of grid points used for the Zernike polynomilas

    :return:
        Z: (np.ndarray(ndim=3, dtype=np.float64)): Cube of Zernike modes
    
    :typical call:
        calc_zcube(n_zpol = self.config.p_z_ab.num_zpol, 
                   grid_p = self.config.p_geom._mpupil.shape[0])
        
    Update of function on Tue Feb  5 10:20:22 2019
    
    @author: hadorn
    
    - Every Zernicke polynomial can now be calculated 
        
    """
    # check if gird_p is a positive int
    if (type(grid_p) != int) and (type(grid_p) != np.int64):
        raise TypeError("grid_p must be a POSITIVE integer.")
    if grid_p <= 0:
        raise TypeError("grid_p must be a POSITIVE integer.")
        
   # create Cartesian grid of needed size (number of grid points)
    (x, y) = np.meshgrid(np.linspace(-1, 1, grid_p), np.linspace(-1, 1, grid_p))

    # transform Cartesian to polar grid
    (r, phi) = cart2pol(x, y)

    # allocate memory for evaluation of Zernike polynomials
    Z = np.zeros((n_zpol, grid_p, grid_p))
    
    # define first row of 'cube' of Zernike coeffs
    Z[0,:,:] = np.ones((grid_p, grid_p))
        
    # iteration for computing all Zernike coefficents
    for zpol in range(1,n_zpol+1):
            
        # determine for noll indicies 
        n,m = Noll2NMindex(zpol)
        
        # iteration for Zernike coeffs does only allow n,m as nonnegativ integers. 
        # If m is smaller than zero, the sign of m has to be changed
        if m < 0:
            m = -m 
        
        # initialize coefficient
        R = 0
        
        # iteration for the zpol-th Zernike coefficient  
        for k in range(0,int((n-m)/2)+1):
            
            # sum computationally regulations of Zernike coefficients
            R += pow(-1,k) * binomial(n-k, k) * binomial(n-2*k, (n-m)/2-k) * np.power(r, n-2*k) 
                    
        # determine normalization factor for Zernike coefficient
        N = 0
        
        # calculate normation based on m -> kronecker delta function
        if m == 0:
            N = math.sqrt(2*(n+1)/2)    
        else:
            N = math.sqrt(2*(n+1))
        
        # check if current calculated zernike polynom is even or odd and if m is equals zero 
        # and determine the Zernike Polynomial
        if (zpol % 2) == 0 and zpol != 1:
            
            Z[zpol-1,:,:] = N*R*np.cos(m*phi)
            
        elif m == 0:
            
            Z[zpol-1,:,:] = N*R
        
        elif (zpol % 2) != 0 and zpol != 1:
                   
            Z[zpol-1,:,:] = N*R*np.sin(m*phi)
    
    return Z
    
#%%
def load_variable(f_dir, f_name, mat_ver, var_name):
    """
    Loads user specified variable from mat-file containing the timeseries of
    Zernike coefficients
    
    :parameters:
        f_dir: (string):    directory of mat-file containing Zernike coeff.
        f_name: (string):   file of mate-fil containing Zernike coeff.
        mat_ver:  (string): Version of the mat-file containing Zernike coeff.
        var_name: (string):   Variable-name of variable to be loaded

    :return:
        coeff: (np.ndarray(ndim=2, dypte=float)): timeseries of Zernike
                                   coefficients without time stamps
    
    :typical call:
    load_variable(f_dir = self.config.p_z_ab.file_dir,
                  f_name = self.config.p_z_ab.file_name,
                  mat_ver = self.config.p_z_ab.mat_vers,
                  var_name = self.config.p_z_ab.var_name_c_xxx *OR* var_name_time)
    """
    # load mat-file via scipy.io
    if mat_ver in ["v4", "v6", "v7"]:
        ret = scio.loadmat(f_dir + f_name)
        
    # load mat-file via hdf5storage
    elif mat_ver == "v7.3":
        ret = hdf5storage.loadmat(f_dir + f_name)
        
    # unexpected error
    else:
        raise TypeError("Unexpected Error: mat_vers must be v4 or v6 or v7 or v7.3")
    
    # return requested vector
    return ret[var_name]

#%%
def calc_step(time, it_time, tol = 1e-10):
    """
    Calculates the time steps of time_series 
    
    :parameters:
        time: (np.ndarray[ndim=1, dtype=np.float64): Time stamps (in seconds) 
                                  of timeseries of Zernike coefficients 
        it_time: (float): COMPASS iteration time
        tol: (float64): acceptable numeric tolerance between 
                        * biggest and smallest step in time
                        * integer and its corresponding float representation

    :return:
        step: (float64): Time steps (in seconds) of time_series
    
    :typical call:
    calc_step(time = self.config.p_z_ab.time_series,
              it_time = self.config.p_loop.ittime,
              tol = 1e-10)
    """
    # calc all time steps
    time_diff = np.diff(time)
    
    # check if time steps vary
    if (time_diff.max() - time_diff.min()) > tol:
        raise ArithmeticError("time_series is not equally spaced")
    
    # calc time step
    step = np.mean(time_diff)
    
    # ratio of step and iteration time
    dec = np.float64( step / it_time )
    
    # check if step is smaller than iteration time
    if dec < 1.0:
        raise ArithmeticError("steps of time_series are smaller than the iteration time")
    
    # check if step is not a multiple of iteration time
    if np.abs( dec - np.round(dec) ) > tol:
        raise ArithmeticError("steps of time_series are not a multiple of the iteration time")
    
    # return step
    return step
 
#%%
def cart2pol(x, y):
    """
    Transforms corresponding elements of the two-dimensional Cartesian 
    coordinate arrays x and y into polar coordinates r and phi
    
    :parameters:
        x: (np.ndarray(ndim=2, dypte=float)): Cartesian coordinate array
        y: (np.ndarray(ndim=2, dypte=float)): Cartesian coordinate array
        example: (x, y) = np.meshgrid(...)

    :return:
        r: (np.ndarray(ndim=2, dypte=float)): polar coordinate array
        phi: (np.ndarray(ndim=2, dypte=float)): polar coordinate array
    """
    # check if r and phi are 2D-matrices
    if len(x.shape) != 2:
        raise TypeError("r must be a ndarray with exactly 2 dimensions")
    if len(y.shape) != 2:
        raise TypeError("phi must be a ndarray with exactly 2 dimensions")
    
    # transformation
    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    
    # return polar coordinates
    return(r, phi)

#%%
def pol2cart(r, phi):
    """
    Transforms corresponding elements of the polar coordinate arrays r and phi 
    to two-dimensional Cartesian coordinates (x and y)
    
    :parameters:
        r: (np.ndarray(ndim=2, dypte=float)): polar coordinate array
        phi: (np.ndarray(ndim=2, dypte=float)): polar coordinate array

    :return:
        x: (np.ndarray(ndim=2, dypte=float)): Cartesian coordinate array
        y: (np.ndarray(ndim=2, dypte=float)): Cartesian coordinate array
        
    """
    # check if r and phi are 2D-matrices
    if len(r.shape) != 2:
        raise TypeError("r must be a ndarray with exactly 2 dimensions")
    if len(phi.shape) != 2:
        raise TypeError("phi must be a ndarray with exactly 2 dimensions")
    
    # transformation
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    
    # return Cartesian coordinates
    return(x, y)

#%%    
def calc_phase_screen(cube, coeff, n_zpol):
    """
    Calcualtes the pupil phase screen at a given time step. Uses the Zernike
    coefficients of respective time step.
    
    :parameters:
        cuebe: (np.ndarray[ndim=1, dtype=np.float64): Time stamps (in seconds) 
                                  of timeseries of Zernike coefficients 
        coeff: (np.ndarray(ndim=1, dypte=np.float64)): timeslice of Zernike
                                   coefficients without time stamps
        n_zpol: (integer):  Number of Zernike polynomials

    :return:
        phase: (np.ndarray[ndim=1, dtype=np.float32): 
    
    :typical call:
        calc_phase_screen(cube = self.config.p_z_ab.zcube_spup *OR* zcube_mpup,
                          coeff = self.config.p_z_ab.coeff_xxx[0,:],
                          n_zpol = self.config.p_z_ab.num_zpol)
    """
    # check if number of coefficients and Zernike polynomials is equal
    if len(coeff) != len(cube[:,0,0]):
        raise ArithmeticError("The number of coefficients and Zernike polynomials must be equal.")
    
    # chekc if number of coefficients and num_zpol is equal
    if n_zpol != len(coeff):
        raise ArithmeticError("The number of coefficients, Zernike polynomials and " + 
                              "num_zpol must be equal.")
    
    # allocate phase-screen-matrix
    phase = np.zeros(shape = cube[0,:,:].shape, dtype = np.float64)
    
    # add every scaled Zernike polynomial to the phase screen matrix
    for k in range(len(coeff)):
        phase = phase + coeff[k] * cube[k,:,:]
    
    # return phase screen (32-bit float)
    return np.float32(phase)

#%%
def calc_trunc(spupil_size, tel_diam, diam):
    """
    Calcualtes the number of pixel a specific diameter would have in the pupil.
    
    :parameters:
        spupil_size: (integer): pupil diameter in pixels
        tel_diam: (float): pupil diameter in meter
        diam: (np.float64): user specified diameter in meter

    :return:
        pix_diam: (np.int64): number of pixel a specific diameter in pupil
    
    :typical call:
        z_aber.calc_trunc(spupil_size = self.config.p_geom.pupdiam,
                          tel_diam = self.config.p_z_ab.pup_diam,
                          diam = self.config.p_z_ab.diam_data)
    """
    # calc pupil scale (m per pix)
    scale = np.float64(tel_diam / spupil_size)
    
    # convert diam (meter) to pixel in pupil
    pix_diam = diam / scale
    
    # round pix_diam to next even int
    pix_diam = np.ceil( np.floor(pix_diam) / 2. ) * 2.
    
    # return int
    return np.int64(pix_diam)
    
#%%
def init_z_aber(sim):
    """ Initializes the custom zernike aberrations (in pupil)
    
    parameters:
        sim: (Simulator) : simulator object (in simulator.py = self)
    """
    
    if sim.config.p_z_ab.include_zab:
        
        # set pup_diam if wildcard is used
        if sim.config.p_z_ab.pup_diam == -1.0:
            # pup_diam = diameter of telescope pupil
            sim.config.p_z_ab.set_pup_diam(sim.config.p_tel.diam)
        
        # if pup_diam is given in meter
        if sim.config.p_z_ab.pup_diam > 0.0:
            # raise error if pup_diam >= diam_data
            if sim.config.p_z_ab.pup_diam >= sim.config.p_z_ab.diam_data:
                raise ArithmeticError("The telescope diameter (pup_diam) must" +
                                      "be smaller than the data diameter (diam_data)") 
            
            # calculate diameter of data (in pixel)
            diam_data_pix = calc_trunc(
                    spupil_size = sim.config.p_geom.pupdiam,
                    tel_diam = sim.config.p_z_ab.pup_diam,
                    diam = sim.config.p_z_ab.diam_data)
            
            # calculate cube of untruncated Zernike polynomials (size of data)
            zcube_untrunc = calc_zcube(
                    n_zpol = sim.config.p_z_ab.num_zpol, 
                    grid_p = diam_data_pix)
            
            # aberration is LARGER than mpupil
            if diam_data_pix > sim.config.p_geom._n:
                # Calculate matrix-padding on one side (mpupil and spupil)
                pix_diff_mp = np.int((diam_data_pix - sim.config.p_geom._n) / 2)
                pix_diff_sp = np.int((diam_data_pix - sim.config.p_geom.pupdiam) / 2)
                
                # truncate cube of Zernike modes for spupil (zcube_spup)
                sim.config.p_z_ab.set_zcube_spup( 
                        zcube_untrunc[:, pix_diff_sp:-pix_diff_sp, pix_diff_sp:-pix_diff_sp] )
        
                # truncate cube of Zernike modes for mpupil (zcube_mpup)
                sim.config.p_z_ab.set_zcube_mpup( 
                        zcube_untrunc[:, pix_diff_mp:-pix_diff_mp, pix_diff_mp:-pix_diff_mp] )
            
            # aberration is AS LARGE AS mpupil
            elif diam_data_pix == sim.config.p_geom._n:
                # Calculate matrix-padding on one side (mpupil and spupil)
                pix_diff_sp = np.int((diam_data_pix - sim.config.p_geom.pupdiam) / 2)
                
                # truncate cube of Zernike modes for spupil (zcube_spup)
                sim.config.p_z_ab.set_zcube_spup( 
                        zcube_untrunc[:, pix_diff_sp:-pix_diff_sp, pix_diff_sp:-pix_diff_sp] )
        
                # cube of Zernike modes for mpupil (zcube_mpup)
                sim.config.p_z_ab.set_zcube_mpup( zcube_untrunc )
            
            # aberration is SMALLER than mpupil
            else:
                raise ArithmeticError("The diameter of the Zernike aberration " + 
                                      "must be as large as the mpupil.")
            
        
        # setup if wildcard is used
        elif sim.config.p_z_ab.pup_diam == -2.0:
            # calculate cube of untruncated Zernike polynomials (size of data)
            zcube_untrunc = calc_zcube(
                    n_zpol = sim.config.p_z_ab.num_zpol, 
                    grid_p = sim.config.p_geom._n)
            
            # Calculate matrix-padding on one side (spupil)
            pix_diff_sp = np.int((sim.config.p_geom._n - sim.config.p_geom.pupdiam) / 2)
            
            # truncate cube of Zernike modes for spupil (zcube_spup)
            sim.config.p_z_ab.set_zcube_spup( 
                    zcube_untrunc[:, pix_diff_sp:-pix_diff_sp, pix_diff_sp:-pix_diff_sp] )
    
            # cube of Zernike modes for mpupil (zcube_mpup)
            sim.config.p_z_ab.set_zcube_mpup( zcube_untrunc )
        
        # error if sim.config.p_z_ab.pup_diam is not > 0 or = -2
        else:
            raise TypeError("pup_diam must be a float >0. or =-1. or =-2.")
        
        # load timeseries of Zernike coefficients (coeff_wfs and coeff_sci) from mat-file
        sim.config.p_z_ab.set_coeff_wfs( load_variable(
                f_dir = sim.config.p_z_ab.file_dir,
                f_name = sim.config.p_z_ab.file_name,
                mat_ver = sim.config.p_z_ab.mat_vers,
                var_name = sim.config.p_z_ab.var_name_c_wfs) )
        
        sim.config.p_z_ab.set_coeff_sci( load_variable(
                f_dir = sim.config.p_z_ab.file_dir,
                f_name = sim.config.p_z_ab.file_name,
                mat_ver = sim.config.p_z_ab.mat_vers,
                var_name = sim.config.p_z_ab.var_name_c_sci) )
        
        # load time stamps coefficients (time_series) from mat-file
        sim.config.p_z_ab.set_time_series( load_variable(
                f_dir = sim.config.p_z_ab.file_dir,
                f_name = sim.config.p_z_ab.file_name,
                mat_ver = sim.config.p_z_ab.mat_vers,
                var_name = sim.config.p_z_ab.var_name_time).ravel() )
        
        # calculate time steps of time_series (step)
        sim.config.p_z_ab.set_step( calc_step(
                time = sim.config.p_z_ab.time_series, 
                it_time = sim.config.p_loop.ittime,
                tol = 1e-12) )
        
        # calculate decimation of iteration time to steps of time_series
        sim.config.p_z_ab.set_dec( int( 
                np.round(sim.config.p_z_ab.step / sim.config.p_loop.ittime) ) )
        
        # calculate and set phase-screen for spupil
        if (sim.config.p_z_ab.include_path == 1) or (sim.config.p_z_ab.include_path == 3):
            sim.tel.set_phase_ab_M1( calc_phase_screen(
                    cube = sim.config.p_z_ab.zcube_spup,
                    coeff = sim.config.p_z_ab.coeff_sci[0,:],
                    n_zpol = sim.config.p_z_ab.num_zpol) )
        
        # calculate and set phase-screen for mpupil
        if (sim.config.p_z_ab.include_path == 2) or (sim.config.p_z_ab.include_path == 3):
            sim.tel.set_phase_ab_M1_m( calc_phase_screen(
                    cube = sim.config.p_z_ab.zcube_mpup,
                    coeff = sim.config.p_z_ab.coeff_wfs[0,:],
                    n_zpol = sim.config.p_z_ab.num_zpol) ) 
        
        # print status of aberrations
        print()
        print("*-------------------------------")
        print("CUSTOM ZERNIKE ABERRATIONS")
        print("status: enabled")
        print("source file: %s" % sim.config.p_z_ab.file_dir + sim.config.p_z_ab.file_name)
        print("number of modes: " + str(sim.config.p_z_ab.num_zpol))
        
        if sim.config.p_z_ab.include_path == 0:
            print("inclusion: not included")
        elif sim.config.p_z_ab.include_path == 1:
            print("inclusion: science (target) path")
        elif sim.config.p_z_ab.include_path == 2:
            print("inclusion: analytic (WFS) path")
        elif sim.config.p_z_ab.include_path == 3:
            print("inclusion: science (target) and analytic (WFS) path")
        else:
            raise TypeError("Something strange happened. " + 
                            "include_path must be 0, 1, 2 or 3.")
        
        if sim.config.p_z_ab.pup_diam > 0.0:
            print("telescope diameter: " + str(sim.config.p_z_ab.pup_diam))
        elif sim.config.p_z_ab.pup_diam == -2.0:
            print("telescope diameter: data fitted to simulation")
        else:
            raise TypeError("Oops, something strange happened. " + 
                            "pup_diam must be > 0.0 or -1.0 or -2.0.")
        print("*-------------------------------")
    
    # nothing to initialize
    else:
        print()
        print("*-------------------------------")
        print("CUSTOM ZERNIKE ABERRATIONS")
        print("status: disabled")
        print("*-------------------------------")
        
#%%
def binomial(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke.
    See http://stackoverflow.com/questions/3025162/statistics-combinations-in-python
    
    parameters:
        n: (integer) : binomial coefficient
        k: (integer) : binomial coefficient
    """
    
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        
        # convert k and n to integer
        n = int(n)
        k = int(k) 
        
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

#%%
def Noll2NMindex(j):
    """
    Converts Noll index to n-m-index
    
    parameters:
        j: (integer) : Noll's index)
    """
    
    n = np.floor( np.sqrt(2*j - 1) + 0.5 ) - 1;

    # calcualte absolute value of azimuthal index
    if np.mod(n,2) == 0:
        # even radial index
        m = 2 * np.floor( (2*j + 1 - n*(n+1)) / 4 );
    else:
        # odd radial index
        m = 2 * np.floor( (2*(j + 1) - n*(n+1)) / 4 ) - 1 ;


    # set sign of azimuthal index
    if np.mod(j,2) != 0:
        m = -m;

    return n,m        
