import pytest
import os
import subprocess

def test_grl40_bm5_paleo17a_BH0_beta1_scalar():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_BH0 -iv beta1 -jsf inputs.json -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_BH0_beta2_scalar():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_BH0 -iv beta2 -jsf inputs.json -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_BH0_Pmax_scalar():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_BH0 -iv Pmax -jsf inputs.json -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_BH0_mu_scalar():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_BH0 -iv mu -jsf inputs.json -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_BH0_s_stat_scalar():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_BH0 -iv s_stat -jsf inputs.json -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_BH0_gamma_s():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_BH0 -iv gamma_s -delta 1.e-4 -jsf inputs.json', shell = True, check = True)

def test_grl40_bm5_paleo17a_BH0_BM5_gamma_s():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_BH0_BM5 -iv gamma_s -delta 1.e-4 -jsf inputs.json', shell = True, check = True)

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

def test_grl40_bm5_paleo17a_BH0_H():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_BH0 -jsf inputs.json', shell = True, check = True)

def test_grl40_bm5_paleo17a_BH0_BM5_H():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_BH0_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_BH0_BM5 -jsf inputs.json', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_H():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -jsf inputs.json', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_BM5_H():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_grl16_bm5_ss25ka_BM5 -jsf inputs.json', shell = True, check = True)

def test_repo_ant40_b2_ss25ka_q_geo():
	subprocess.run (['cp', 'headers/sico_specs_repo_ant40_b2_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_ant40_b2_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_ant40_b2_ss25ka -dom ant -iv q_geo -jsf inputs.json', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_q_geo():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -iv q_geo -jsf inputs.json', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_BM5_q_geo():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_grl16_bm5_ss25ka_BM5 -iv q_geo -jsf inputs.json', shell = True, check = True)

def test_repo_ant40_b2_ss25ka_temp_c():
	subprocess.run (['cp', 'headers/sico_specs_repo_ant40_b2_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_ant40_b2_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_ant40_b2_ss25ka -dom ant -iv temp_c -dim 3 -z 40 -delta 5.e-2 -jsf inputs.json', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_temp_c():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -iv temp_c -dim 3 -z 40 -jsf inputs.json', shell = True, check = True)