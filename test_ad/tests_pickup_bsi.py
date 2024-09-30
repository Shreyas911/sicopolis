import pytest
import os
import subprocess

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup -iv c_slide_init -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup -iv delta_tda_const -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup -iv beta1 -delta 3.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup -iv beta2 -delta 7.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup -iv Pmax -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup -iv mu -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup -iv s_stat -delta 5.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup_gamma_s_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus1ka_pickup -iv gamma_s -delta 7.e-3 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup -iv c_slide_init -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup -iv delta_tda_const -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup -iv beta1 -delta 3.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup -iv beta2 -delta 7.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup -iv Pmax -delta 1.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup -iv mu -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup -iv s_stat -delta 5.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup_gamma_s_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus1ka_pickup -iv gamma_s -delta 7.e-3 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup -iv c_slide_init -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup -iv delta_tda_const -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup -iv beta1 -delta 3.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup -iv beta2 -delta 7.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup -iv Pmax -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup -iv mu -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup -iv s_stat -delta 5.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup_gamma_s_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus1ka_pickup -iv gamma_s -delta 7.e-3 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup -iv c_slide_init -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup -iv delta_tda_const -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup -iv beta1 -delta 3.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup -iv beta2 -delta 7.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup -iv Pmax -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup -iv mu -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup -iv s_stat -delta 5.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup_gamma_s_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus1ka_pickup -iv gamma_s -delta 7.e-3 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup -iv c_slide_init -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup -iv delta_tda_const -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup -iv beta1 -delta 3.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup -iv beta2 -delta 7.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup -iv Pmax -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup -iv mu -delta 5.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup -iv s_stat -delta 5.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup_gamma_s_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_minus5ka_pickup -iv gamma_s -delta 7.e-3 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup -iv c_slide_init -delta 5.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup -iv delta_tda_const -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup -iv beta1 -delta 3.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup -iv beta2 -delta 7.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup -iv Pmax -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup -iv mu -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup -iv s_stat -delta 5.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup_gamma_s_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_minus5ka_pickup -iv gamma_s -delta 7.e-3 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup -iv c_slide_init -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup -iv delta_tda_const -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup -iv beta1 -delta 3.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup -iv beta2 -delta 7.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup -iv Pmax -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup -iv mu -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup -iv s_stat -delta 5.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup_gamma_s_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_minus5ka_pickup -iv gamma_s -delta 7.e-3 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup -iv c_slide_init -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup -iv delta_tda_const -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup -iv beta1 -delta 3.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup -iv beta2 -delta 7.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup -iv Pmax -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup -iv mu -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup -iv s_stat -delta 5.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup_gamma_s_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_minus5ka_pickup -iv gamma_s -delta 7.e-3 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

