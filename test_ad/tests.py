import pytest
import os
import subprocess

def test_repo_ant64_bm3_ss25ka_H():
	subprocess.run (['cp', 'headers/sico_specs_repo_ant64_bm3_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_ant64_bm3_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_ant64_bm3_ss25ka -dom ant -jsf inputs.json', shell = True, check = True)

def test_repo_ant40_b2_ss25ka_H():
	subprocess.run (['cp', 'headers/sico_specs_repo_ant40_b2_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_ant40_b2_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_ant40_b2_ss25ka -dom ant -jsf inputs.json', shell = True, check = True)

def test_repo_ant64_b2_future09_ctrl_H():
	subprocess.run (['cp', 'headers/sico_specs_repo_ant64_b2_future09_ctrl.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_ant64_b2_future09_ctrl.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_ant64_b2_future09_ctrl -dom ant -delta 5.e-2 -jsf inputs.json', shell = True, check = True)

def test_grl10_bm5_paleo17a_BH0_H():
	subprocess.run (['cp', 'headers/sico_specs_grl10_bm5_paleo17a_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl10_bm5_paleo17a_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl10_bm5_paleo17a_BH0 -jsf inputs.json', shell = True, check = True)

def test_grl40_bm5_paleo17a_BH0_H():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_BH0 -jsf inputs.json', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_H():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -jsf inputs.json', shell = True, check = True)

def test_repo_ant40_b2_ss25ka_q_geo():
	subprocess.run (['cp', 'headers/sico_specs_repo_ant40_b2_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_ant40_b2_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_ant40_b2_ss25ka -dom ant -iv q_geo -jsf inputs.json', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_q_geo():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -iv q_geo -jsf inputs.json', shell = True, check = True)

def test_repo_ant40_b2_ss25ka_temp_c():
	subprocess.run (['cp', 'headers/sico_specs_repo_ant40_b2_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_ant40_b2_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_ant40_b2_ss25ka -dom ant -iv temp_c -dim 3 -z 40 -delta 5.e-2 -jsf inputs.json', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_temp_c():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -iv temp_c -dim 3 -z 40 -jsf inputs.json', shell = True, check = True)