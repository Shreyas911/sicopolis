import pytest
import os
import subprocess

def test_v5_ant40_ss25ka():
	subprocess.run('python tapenade_config.py -head v5_ant40_ss25ka -dom ant -dv fc -iv H -delta 1.e-3 --travis', capture_output = True, shell = True, check = True)

def test_v5_ant64_b2_future09_ctrl():
	subprocess.run('python tapenade_config.py -head v5_ant64_b2_future09_ctrl -dom ant -dv fc -iv H -delta 1.e-3 --travis', capture_output = True, shell = True, check = True)

def test_v5_grl20_ss25ka():	
	subprocess.run('python tapenade_config.py -head v5_grl20_ss25ka -dom grl -dv fc -iv H -delta 1.e-3 --travis', capture_output = True, shell = True, check = True)
