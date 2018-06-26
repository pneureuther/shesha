"""
Mat Logger Utility

Data logging via mat-file

@author: neureuther
"""

import sys
import os
import hdf5storage


def mat_file_init(log_mat, mat_file_name, mat_file_dir, mat_var, bool_log_init):
    """
    Initializes the logging of variables in a mat-file
    
    :parameters:
        log_mat: (bool) : variables are saved to mat-file if true
            = self.config.p_loop.log_mat
        mat_file_name: (string) : file-name of mat-file
            = self.config.p_loop.mat_file_name
        mat_file_dir: (string) : directory of mat-file
            = self.config.p_loop.mat_file_dir
        mat_var: (list of strings) : list of variable names used in mat-file
            = self.config.p_loop.mat_var
        bool_log_init: (bool) : variables after initialization are saved
            = self.config.p_loop.log_init
    """
    
    # Initialize File
    if log_mat:
        print()
        print("*-------------------------------")
        print("DATA-LOGGING")
        print("status: enabled")
        if bool_log_init:
            print("log init-values: enabled")
        else:
            print("log init-values: disabled")
        
        # create complete file-name
        mat_file = mat_file_dir + mat_file_name
        
        # check if mat-file already exists
        if os.path.isfile(mat_file):
            print("file name: %s" % mat_file_name)
            print("file location: %s" % mat_file)
            reply = str(input("This file already exists. Do you want to overwrite it? (y/n): ")).lower().strip()
            while True:
                if reply[:1] == 'y':
                    os.remove(mat_file)
                    print("%s overwritten" % mat_file_name)
                    break
                elif reply[:1] == 'n':
                    raise ValueError("%s already exists and was NOT overwritten"
                                     % mat_file_name)
                    break
                else:
                    reply = str(input("Invalid answer. Do you want to overwrite it? (y/n): ")).lower().strip()
        else:
            print("%s created" % mat_file_name)
            print("file location: %s" % mat_file)
        
        # create EMPTY mat-file
        hdf5storage.savemat(mat_file, {}, appendmat = False, format = '7.3',
                            action_for_matlab_incompatible = 'error')
        
        # printing names of logged variables
        print("logging variables:")
        for k in range(len(mat_var)):
            print("   * %s" % mat_var[k])
        print("*-------------------------------")
        print()
    
    # Nothing to initialize
    else:
        print()
        print("*-------------------------------")
        print("DATA-LOGGING")
        print("status: disabled")
        print("*-------------------------------")
        print()
        
def mat_file_write(log_mat, mat_file_name, mat_file_dir, mat_num_dig,
                   mat_var, mat_cmds, index, sim):
    """
    Writes new data to the specified mat-file
    
    :parameters:
        log_mat: (bool) : variables are saved to mat-file if true
            = self.config.p_loop.log_mat
        mat_file_name: (string) : file-name of mat-file
            = self.config.p_loop.mat_file_name
        mat_file_dir: (string) : directory of mat-file
            = self.config.p_loop.mat_file_dir
        mat_num_dig: (int) : number of digits used in the mat-file variables
            = self.config.p_loop.mat_num_dig
        mat_var: (list of strings) : list of variable names used in mat-file
            = self.config.p_loop.mat_var
        mat_cmds: (list of strings) : list of COMPASS commands used to retrievedesiered data
            = self.config.p_loop.mat_cmds
        index: (int) : iteration index of optical calculation
    """
    
    # Creating a list of strings:
    # ["Alpha001", "Bravo001", etc.] (mat_var with appended index)
    num_setting = '{:0%d}' % mat_num_dig
    list_names = [None]*len(mat_var)

    for k in range(len(mat_var)):
        list_names[k] = mat_var[k] + num_setting.format(index+1)
    
    # Creating a list of variables to be saved in mat-file
    list_vars = [None]*len(mat_cmds)
    
    for k in range(len(mat_cmds)):
        cmd = "list_vars[%i]" % k
        cmd = cmd + " = " + mat_cmds[k]
        exec(cmd)
    
    # Create dictionary with variables to be saved
    mat_dict = dict(zip(list_names, list_vars))
    
    # Save desired variables to mat-file
    hdf5storage.writes(mat_dict, filename = mat_file_dir + mat_file_name,
                       matlab_compatible = True)
    
    