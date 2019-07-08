#!/usr/bin/env python
"""script test to simulate a closed loop

Usage:
  closed_loop.py <parameters_filename>

with 'parameters_filename' the path to the parameters file
"""

from docopt import docopt

if __name__ == "__main__":
    import shesha.sim

    arguments = docopt(__doc__)
    param_file = arguments["<parameters_filename>"]

    from shesha.sim.simulatorCE import SimulatorCE

    sim = SimulatorCE(param_file)#, True)

    sim.init_sim()
    sim.loop(n = sim.config.p_loop.niter, save_index = 0, phase_index = 1)#, do_control = False, see_atmos = False)
