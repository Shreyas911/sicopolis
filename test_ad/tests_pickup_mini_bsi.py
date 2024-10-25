import pytest
import os
import subprocess

def test_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini_temp_c_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini -iv temp_c -dim 3 -z 40 -low 0.01 -jsf inputs.json -asi 0 -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini_age_c_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini -iv age_c -dim 3 -z 40 -jsf inputs.json -asi 0 -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini_H_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini -iv H -delta 1.e-1 -low 0.01 -jsf inputs_pickup.json -asi 0 -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini_q_geo_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini -iv q_geo -jsf inputs_pickup.json -asi 0 -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini -iv c_slide_init -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini -iv delta_tda_const -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 0', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini -iv beta1 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini -iv beta2 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini -iv Pmax -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini -iv mu -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini -iv s_stat -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini_gamma_s_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini -iv gamma_s -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini_c_dis_da_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_AC_m1ka_pkp_mini -iv c_dis_da -delta 1.6e+2 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini_temp_c_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini -iv temp_c -dim 3 -z 40 -low 0.01 -jsf inputs.json -asi 0 -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini_age_c_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini -iv age_c -dim 3 -z 40 -jsf inputs.json -asi 0 -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini_H_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini -iv H -delta 1.e-1 -low 0.01 -jsf inputs_pickup.json -asi 0 -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini_q_geo_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini -iv q_geo -jsf inputs_pickup.json -asi 0 -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini -iv c_slide_init -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini -iv delta_tda_const -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 0', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini -iv beta1 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini -iv beta2 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini -iv Pmax -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini -iv mu -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini -iv s_stat -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini_gamma_s_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini -iv gamma_s -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini_c_dis_da_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_BM5_m1ka_pkp_mini -iv c_dis_da -delta 1.6e+2 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini_temp_c_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini -iv temp_c -dim 3 -z 40 -low 0.01 -jsf inputs.json -asi 0 -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini_age_c_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini -iv age_c -dim 3 -z 40 -jsf inputs.json -asi 0 -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini_H_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini -iv H -delta 1.e-1 -low 0.01 -jsf inputs_pickup.json -asi 0 -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini_q_geo_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini -iv q_geo -jsf inputs_pickup.json -asi 0 -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini -iv c_slide_init -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini -iv delta_tda_const -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 0', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini -iv beta1 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini -iv beta2 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini -iv Pmax -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini -iv mu -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini -iv s_stat -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini_gamma_s_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini -iv gamma_s -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini_c_dis_da_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FAC_m1ka_pkp_mini -iv c_dis_da -delta 1.6e+2 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini_temp_c_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini -iv temp_c -dim 3 -z 40 -low 0.01 -jsf inputs.json -asi 0 -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini_age_c_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini -iv age_c -dim 3 -z 40 -jsf inputs.json -asi 0 -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini_H_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini -iv H -delta 1.e-1 -low 0.01 -jsf inputs_pickup.json -asi 0 -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini_q_geo_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini -iv q_geo -jsf inputs_pickup.json -asi 0 -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini -iv c_slide_init -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini -iv delta_tda_const -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 0', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini -iv beta1 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini -iv beta2 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini -iv Pmax -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini -iv mu -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini -iv s_stat -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini_gamma_s_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini -iv gamma_s -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini_c_dis_da_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_FBM5_m1ka_pkp_mini -iv c_dis_da -delta 1.6e+2 -jsf inputs_pickup.json -asi 0 -lbfs scalar -bia 1', shell = True, check = True)

