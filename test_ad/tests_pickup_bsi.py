import pytest
import os
import subprocess

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp -iv c_slide_init -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp -iv delta_tda_const -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp -iv beta1 -delta 3.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp -iv beta2 -delta 7.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp -iv Pmax -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp -iv mu -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp -iv s_stat -delta 5.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp_gamma_s_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FAC_m1ka_pkp -iv gamma_s -delta 7.e-3 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp -iv c_slide_init -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp -iv delta_tda_const -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp -iv beta1 -delta 3.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp -iv beta2 -delta 7.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp -iv Pmax -delta 1.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp -iv mu -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp -iv s_stat -delta 5.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp_gamma_s_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_AC_m1ka_pkp -iv gamma_s -delta 7.e-3 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp -iv c_slide_init -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp -iv delta_tda_const -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp -iv beta1 -delta 3.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp -iv beta2 -delta 7.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp -iv Pmax -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp -iv mu -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp -iv s_stat -delta 5.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp_gamma_s_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_FBM5_m1ka_pkp -iv gamma_s -delta 7.e-3 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp_c_slide_init_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp -iv c_slide_init -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp_delta_tda_const_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp -iv delta_tda_const -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp_beta1_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp -iv beta1 -delta 3.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp_beta2_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp -iv beta2 -delta 7.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp_Pmax_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp -iv Pmax -delta 5.e-2 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp_mu_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp -iv mu -delta 1.e+0 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp_s_stat_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp -iv s_stat -delta 5.e-1 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

def test_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp_gamma_s_scalar_bsi():
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp.h', '../headers'])
	subprocess.run (['cp', 'headers/sico_specs_grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp.h', '../src/sico_specs.h'])
	subprocess.run('python3 tapenade_config.py -head grl40_bm5_paleo17a_CT4_BH0_tunedCS_BM5_m1ka_pkp -iv gamma_s -delta 7.e-3 -jsf inputs_pickup.json -asi 0 -lbfs scalar', shell = True, check = True)

