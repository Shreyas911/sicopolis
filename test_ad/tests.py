import pytest
import os
import subprocess

def test_v5_ant40_ss25ka_H():
	subprocess.run (['cp', '../runs/headers/sico_specs_v5_ant40_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head v5_ant40_ss25ka -dom ant -dv fc -iv H -delta 1.e-3 --travis', shell = True, check = True)

def test_v5_ant64_b2_future09_ctrl_H():
	subprocess.run (['cp', '../runs/headers/sico_specs_v5_ant64_b2_future09_ctrl.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head v5_ant64_b2_future09_ctrl -dom ant -dv fc -iv H -delta 5.e-2 --travis', shell = True, check = True)

def test_v5_grl20_ss25ka_H():	
	subprocess.run (['cp', '../runs/headers/sico_specs_v5_grl20_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head v5_grl20_ss25ka -dom grl -dv fc -iv H -delta 1.e-3 --travis', shell = True, check = True)

def test_v5_ant40_ss25ka_q_geo():
	subprocess.run (['cp', '../runs/headers/sico_specs_v5_ant40_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head v5_ant40_ss25ka -dom ant -dv fc -iv q_geo -delta 1.e-3 --travis', shell = True, check = True)

def test_v5_grl20_ss25ka_q_geo():
        subprocess.run (['cp', '../runs/headers/sico_specs_v5_grl20_ss25ka.h', '../src/sico_specs.h'])
        subprocess.run('python3 tapenade_config.py -head v5_grl20_ss25ka -dom grl -dv fc -iv q_geo -delta 1.e-3 --travis', shell = True, check = True)
