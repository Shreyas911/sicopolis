import pytest
import subprocess

@pytest.mark.parametrize("header_filename, head, iv, delta, low, dim, z, jsf, asi, lbfs, bia, dom", [
    # delta_tda_const cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "delta_tda_const", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "delta_tda_const", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_BM5.h", "repo_grl16_bm5_ss25ka_BM5", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_ant64_b2_future09_ctrl.h", "repo_ant64_b2_future09_ctrl", "delta_tda_const", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, "ant"),
    ("sico_specs_repo_ant64_bm3_ss25ka.h", "repo_ant64_bm3_ss25ka", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, "ant"),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "delta_tda_const", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "delta_tda_const", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_SVC.h", "repo_grl16_bm5_ss25ka_CM3_SVC", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),

    # beta1 cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "beta1", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "beta1", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "beta1", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "beta1", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "beta1", "1.e-1", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "beta1", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "beta1", "1.e-1", None, None, None, "inputs.json", 1, "scalar", 1, None),

    # beta2 cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),

    # Pmax cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "Pmax", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "Pmax", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "Pmax", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "Pmax", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "Pmax", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "Pmax", None, None, None, None, "inputs.json", 1, "scalar", 1, None),

    # mu cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "mu", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "mu", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "mu", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "mu", "1.e-1", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "mu", "1.e-1", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "mu", "1.e-1", None, None, None, "inputs.json", 1, "scalar", 1, None),

    # gamma_s cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "gamma_s", "1.e-4", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "gamma_s", "1.e-4", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "gamma_s", "1.e-4", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "gamma_s", "1.e-4", "0.01", None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "gamma_s", "1.e-4", "0.01", None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "gamma_s", "1.e-4", "0.01", None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "gamma_s", "1.e-4", "0.01", None, None, "inputs.json", 1, "scalar", 1, None),

    # H cases
    ("sico_specs_repo_ant64_bm3_ss25ka.h", "repo_ant64_bm3_ss25ka", "H", None, None, None, None, "inputs.json", 1, "scalar", 1, "ant"),
    ("sico_specs_repo_ant64_b2_future09_ctrl.h", "repo_ant64_b2_future09_ctrl", "H", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, "ant"),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "H", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "H", "1.e-4", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "H", "1.e-4", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "H", "1.e-4", "0.01", None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "H", "1.e-4", "0.01", None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "H", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_BM5.h", "repo_grl16_bm5_ss25ka_BM5", "H", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "H", None, "0.01", None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "H", "1.e-4", "0.01", None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "H", "1.e-4", "0.01", None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_SVC.h", "repo_grl16_bm5_ss25ka_CM3_SVC", "H", None, "0.01", None, None, "inputs.json", 1, "scalar", 1, None),

    # q_geo cases
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "q_geo", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_BM5.h", "repo_grl16_bm5_ss25ka_BM5", "q_geo", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "q_geo", None, "0.01", None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_SVC.h", "repo_grl16_bm5_ss25ka_CM3_SVC", "q_geo", None, "0.01", None, None, "inputs.json", 1, "scalar", 1, None),

    # c_slide_init cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "c_slide_init", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_BM5.h", "repo_grl16_bm5_ss25ka_BM5", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "c_slide_init", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "c_slide_init", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_SVC.h", "repo_grl16_bm5_ss25ka_CM3_SVC", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),

    # s_stat cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "s_stat", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "s_stat", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "s_stat", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "s_stat", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "s_stat", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "s_stat", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "s_stat", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),

    # temp_c and age_c cases
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "temp_c", None, None, "3", "40", "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "age_c", None, None, "3", "40", "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "age_c", None, None, "3", "40", "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "age_c", None, None, "3", "40", "inputs.json", 1, "scalar", 1, None),
])
def test_tapenade_config(header_filename, head, iv, delta, low, dim, z, jsf, asi, lbfs, bia, dom):
    subprocess.run(['cp', f'headers/{header_filename}', '../headers'])
    subprocess.run(['cp', f'headers/{header_filename}', '../src/sico_specs.h'])

    cmd = f'python3 tapenade_config.py -head {head} -iv {iv} -jsf {jsf}'

    if delta:
        cmd += f' -delta {delta}'
    if low:
        cmd += f' -low {low}'
    if dim:
        cmd += f' -dim {dim}'
    if z:
        cmd += f' -z {z}'
    if lbfs:
        cmd += f' -lbfs {lbfs}'
    if bia:
        cmd += f' -bia {bia}'
    if asi:
        cmd += f' -asi {asi}'
    if dom:
        cmd += f' -dom {dom}'
    subprocess.run(cmd, shell=True, check=True)
