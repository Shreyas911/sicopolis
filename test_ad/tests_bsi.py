import pytest
import os
import subprocess

def test_grl40_bm5_paleo17a_CT4_BH0_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0 -iv delta_tda_const -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5 -iv delta_tda_const -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5 -iv delta_tda_const -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC -iv delta_tda_const -delta 5.e-2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC -iv delta_tda_const -delta 5.e-2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_grl16_bm5_ss25ka -iv delta_tda_const -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_BM5_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_grl16_bm5_ss25ka_BM5 -iv delta_tda_const -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_CM3_AC_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_grl16_bm5_ss25ka_CM3_AC -iv delta_tda_const -delta 1.e-5 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_repo_ant64_b2_future09_ctrl_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_ant64_b2_future09_ctrl.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_ant64_b2_future09_ctrl.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_ant64_b2_future09_ctrl -iv delta_tda_const -delta 5.e-2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_repo_ant64_bm3_ss25ka_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_ant64_bm3_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_ant64_bm3_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_ant64_bm3_ss25ka -dom ant -iv c_slide_init -delta 1.e-2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5 -iv c_slide_init -delta 1.e-2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5 -iv c_slide_init -delta 1.e-2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -iv c_slide_init -delta 1.e-2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_BM5_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_grl16_bm5_ss25ka_BM5 -iv c_slide_init -delta 1.e-2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_CM3_AC_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_grl16_bm5_ss25ka_CM3_AC -iv c_slide_init -delta 1.e-2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0 -iv beta1 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5 -iv beta1 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC -iv beta1 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5 -iv beta1 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC -iv beta1 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0 -iv beta2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5 -iv beta2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC -iv beta2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5 -iv beta2 -delta 5.e-2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC -iv beta2 -delta 5.e-2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0 -iv Pmax -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5 -iv Pmax -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC -iv Pmax -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5 -iv Pmax -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC -iv Pmax -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0 -iv mu -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5 -iv mu -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC -iv mu -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5 -iv mu -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC -iv mu -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0 -iv s_stat -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5 -iv s_stat -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC -iv s_stat -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5 -iv s_stat -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC -iv s_stat -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_c_dis_da_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0 -iv c_dis_da -delta 1.6e2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_c_dis_da_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5 -iv c_dis_da -delta 1.6e2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_c_dis_da_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC -iv c_dis_da -delta 1.6e2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_c_dis_da_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5 -iv c_dis_da -delta 1.6e2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_c_dis_da_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC -iv c_dis_da -delta 1.6e2 -jsf inputs.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_gamma_s_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0 -iv gamma_s -delta 1.e-4 -jsf inputs.json -asi 0', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_gamma_s_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5 -iv gamma_s -delta 1.e-4 -jsf inputs.json -asi 0', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_gamma_s_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5 -iv gamma_s -delta 1.e-4 -low 0.05 -jsf inputs.json -asi 0', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_gamma_s_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC -iv gamma_s -delta 1.e-4 -low 0.01 -jsf inputs.json -asi 0', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_gamma_s_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC -iv gamma_s -delta 1.e-4 -low 0.01 -jsf inputs.json -asi 0', shell = True, check = True)

def test_repo_ant64_bm3_ss25ka_H_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_ant64_bm3_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_ant64_bm3_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_ant64_bm3_ss25ka -dom ant -jsf inputs.json -asi 0', shell = True, check = True)

def test_repo_ant64_b2_future09_ctrl_H_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_ant64_b2_future09_ctrl.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_ant64_b2_future09_ctrl.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_ant64_b2_future09_ctrl -dom ant -delta 5.e+0 -jsf inputs.json -asi 0', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_H_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0 -delta 1.e-1 -jsf inputs.json -asi 0', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_H_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5 -delta 1.e-4 -jsf inputs.json -asi 0', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_H_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5 -delta 1.e-4 -low 0.01 -jsf inputs.json -asi 0', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_H_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC -delta 1.e-4 -low 0.01 -jsf inputs.json -asi 0', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_H_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC -delta 1.e-4 -low 0.015 -jsf inputs.json -asi 0', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_H_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -low 0.01 -jsf inputs.json -asi 0', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_BM5_H_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_grl16_bm5_ss25ka_BM5 -jsf inputs.json -asi 0', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_CM3_AC_H_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_grl16_bm5_ss25ka_CM3_AC -low 0.01 -jsf inputs.json -asi 0', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_q_geo_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -iv q_geo -jsf inputs.json -asi 0', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_BM5_q_geo_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_BM5.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_BM5.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_grl16_bm5_ss25ka_BM5 -iv q_geo -jsf inputs.json -asi 0', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_CM3_AC_q_geo_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_grl16_bm5_ss25ka_CM3_AC -iv q_geo -delta 1.e-4 -jsf inputs.json -asi 0', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_temp_c_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -iv temp_c -dim 3 -z 40 -low 0.01 -jsf inputs.json -asi 0', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_age_c_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC -iv age_c -dim 3 -z 40 -jsf inputs.json -asi 0', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_age_c_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC -iv age_c -dim 3 -z 40 -jsf inputs.json -asi 0', shell = True, check = True)

def test_repo_grl16_bm5_ss25ka_CM3_AC_age_c_bsi():
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head repo_grl16_bm5_ss25ka_CM3_AC -iv age_c -dim 3 -z 40 -jsf inputs.json -asi 0', shell = True, check = True)

