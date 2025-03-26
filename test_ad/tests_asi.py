import pytest
import subprocess

@pytest.mark.parametrize("header_filename, head, iv, delta, low, dim, z, jsf, asi, lbfs, bia, dom", [

    # RHO_A cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "RHO_A", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "RHO_A", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "RHO_A", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZSC.h", "grl40_bm5_paleo17a_CT4_BH0_ZSC", "RHO_A", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZSC.h", "grl40_bm5_paleo17a_CT4_BH0_FZSC", "RHO_A", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZLC.h", "grl40_bm5_paleo17a_CT4_BH0_ZLC", "RHO_A", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZLC.h", "grl40_bm5_paleo17a_CT4_BH0_FZLC", "RHO_A", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "RHO_A", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "RHO_A", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "RHO_A", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_BM5.h", "repo_grl16_bm5_ss25ka_BM5", "RHO_A", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZSC.h", "repo_grl16_bm5_ss25ka_ZSC", "RHO_A", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZLC.h", "repo_grl16_bm5_ss25ka_ZLC", "RHO_A", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "RHO_A", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_ant64_b2_future09_ctrl.h", "repo_ant64_b2_future09_ctrl", "RHO_A", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, "ant"),
    ("sico_specs_repo_ant64_bm3_ss25ka.h", "repo_ant64_bm3_ss25ka", "RHO_A", None, None, None, None, "inputs.json", 1, "scalar", 1, "ant"),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "RHO_A", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "RHO_A", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),

    # time_lag_asth cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "time_lag_asth", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "time_lag_asth", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "time_lag_asth", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZSC.h", "grl40_bm5_paleo17a_CT4_BH0_ZSC", "time_lag_asth", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZSC.h", "grl40_bm5_paleo17a_CT4_BH0_FZSC", "time_lag_asth", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZLC.h", "grl40_bm5_paleo17a_CT4_BH0_ZLC", "time_lag_asth", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZLC.h", "grl40_bm5_paleo17a_CT4_BH0_FZLC", "time_lag_asth", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "time_lag_asth", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "time_lag_asth", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "time_lag_asth", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_BM5.h", "repo_grl16_bm5_ss25ka_BM5", "time_lag_asth", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZSC.h", "repo_grl16_bm5_ss25ka_ZSC", "time_lag_asth", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZLC.h", "repo_grl16_bm5_ss25ka_ZLC", "time_lag_asth", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "time_lag_asth", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_ant64_b2_future09_ctrl.h", "repo_ant64_b2_future09_ctrl", "time_lag_asth", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, "ant"),
    ("sico_specs_repo_ant64_bm3_ss25ka.h", "repo_ant64_bm3_ss25ka", "time_lag_asth", None, None, None, None, "inputs.json", 1, "scalar", 1, "ant"),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "time_lag_asth", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "time_lag_asth", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),

    # flex_rig_lith cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZSC.h", "grl40_bm5_paleo17a_CT4_BH0_ZSC", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZSC.h", "grl40_bm5_paleo17a_CT4_BH0_FZSC", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZLC.h", "grl40_bm5_paleo17a_CT4_BH0_ZLC", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZLC.h", "grl40_bm5_paleo17a_CT4_BH0_FZLC", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_BM5.h", "repo_grl16_bm5_ss25ka_BM5", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZSC.h", "repo_grl16_bm5_ss25ka_ZSC", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZLC.h", "repo_grl16_bm5_ss25ka_ZLC", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_ant64_b2_future09_ctrl.h", "repo_ant64_b2_future09_ctrl", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, "ant"),
    ("sico_specs_repo_ant64_bm3_ss25ka.h", "repo_ant64_bm3_ss25ka", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, "ant"),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "flex_rig_lith", "1.e-7", None, None, None, "inputs.json", 1, None, 1, None),

    # delta_tda_const cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZSC.h", "grl40_bm5_paleo17a_CT4_BH0_ZSC", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZSC.h", "grl40_bm5_paleo17a_CT4_BH0_FZSC", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZLC.h", "grl40_bm5_paleo17a_CT4_BH0_ZLC", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZLC.h", "grl40_bm5_paleo17a_CT4_BH0_FZLC", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "delta_tda_const", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "delta_tda_const", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_BM5.h", "repo_grl16_bm5_ss25ka_BM5", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZSC.h", "repo_grl16_bm5_ss25ka_ZSC", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZLC.h", "repo_grl16_bm5_ss25ka_ZLC", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_ant64_b2_future09_ctrl.h", "repo_ant64_b2_future09_ctrl", "delta_tda_const", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, "ant"),
    ("sico_specs_repo_ant64_bm3_ss25ka.h", "repo_ant64_bm3_ss25ka", "delta_tda_const", None, None, None, None, "inputs.json", 1, "scalar", 1, "ant"),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "delta_tda_const", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "delta_tda_const", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),

    # c_slide_init cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZSC.h", "grl40_bm5_paleo17a_CT4_BH0_ZSC", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZSC.h", "grl40_bm5_paleo17a_CT4_BH0_FZSC", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZLC.h", "grl40_bm5_paleo17a_CT4_BH0_ZLC", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZLC.h", "grl40_bm5_paleo17a_CT4_BH0_FZLC", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "c_slide_init", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_BM5.h", "repo_grl16_bm5_ss25ka_BM5", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZSC.h", "repo_grl16_bm5_ss25ka_ZSC", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZLC.h", "repo_grl16_bm5_ss25ka_ZLC", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "c_slide_init", "1.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "c_slide_init", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "c_slide_init", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),

    # beta1 cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "beta1", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "beta1", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "beta1", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZSC.h", "grl40_bm5_paleo17a_CT4_BH0_ZSC", "beta1", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZSC.h", "grl40_bm5_paleo17a_CT4_BH0_FZSC", "beta1", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZLC.h", "grl40_bm5_paleo17a_CT4_BH0_ZLC", "beta1", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZLC.h", "grl40_bm5_paleo17a_CT4_BH0_FZLC", "beta1", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "beta1", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "beta1", "1.e-1", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "beta1", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "beta1", "1.e-1", None, None, None, "inputs.json", 1, "scalar", 1, None),

    # beta2 cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZSC.h", "grl40_bm5_paleo17a_CT4_BH0_ZSC", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZSC.h", "grl40_bm5_paleo17a_CT4_BH0_FZSC", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZLC.h", "grl40_bm5_paleo17a_CT4_BH0_ZLC", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZLC.h", "grl40_bm5_paleo17a_CT4_BH0_FZLC", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "beta2", None, None, None, None, "inputs.json", 1, "scalar", 1, None),

    # Pmax cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "Pmax", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "Pmax", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "Pmax", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZSC.h", "grl40_bm5_paleo17a_CT4_BH0_ZSC", "Pmax", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZSC.h", "grl40_bm5_paleo17a_CT4_BH0_FZSC", "Pmax", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZLC.h", "grl40_bm5_paleo17a_CT4_BH0_ZLC", "Pmax", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZLC.h", "grl40_bm5_paleo17a_CT4_BH0_FZLC", "Pmax", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "Pmax", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "Pmax", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "Pmax", None, None, None, None, "inputs.json", 1, "scalar", 1, None),

    # mu cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "mu", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "mu", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "mu", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZSC.h", "grl40_bm5_paleo17a_CT4_BH0_ZSC", "mu", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZSC.h", "grl40_bm5_paleo17a_CT4_BH0_FZSC", "mu", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZLC.h", "grl40_bm5_paleo17a_CT4_BH0_ZLC", "mu", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZLC.h", "grl40_bm5_paleo17a_CT4_BH0_FZLC", "mu", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "mu", "1.e-1", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "mu", "1.e-1", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "mu", "1.e-1", None, None, None, "inputs.json", 1, "scalar", 1, None),

    # s_stat cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "s_stat", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "s_stat", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "s_stat", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZSC.h", "grl40_bm5_paleo17a_CT4_BH0_ZSC", "s_stat", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZSC.h", "grl40_bm5_paleo17a_CT4_BH0_FZSC", "s_stat", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZLC.h", "grl40_bm5_paleo17a_CT4_BH0_ZLC", "s_stat", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZLC.h", "grl40_bm5_paleo17a_CT4_BH0_FZLC", "s_stat", None, None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "s_stat", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "s_stat", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "s_stat", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "s_stat", "5.e-2", None, None, None, "inputs.json", 1, "scalar", 1, None),

    # gamma_s cases
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "gamma_s", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "gamma_s", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "gamma_s", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZSC.h", "grl40_bm5_paleo17a_CT4_BH0_ZSC", "gamma_s", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZSC.h", "grl40_bm5_paleo17a_CT4_BH0_FZSC", "gamma_s", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZLC.h", "grl40_bm5_paleo17a_CT4_BH0_ZLC", "gamma_s", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZLC.h", "grl40_bm5_paleo17a_CT4_BH0_FZLC", "gamma_s", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "gamma_s", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "gamma_s", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "gamma_s", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "gamma_s", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),

    # H cases
    ("sico_specs_repo_ant64_bm3_ss25ka.h", "repo_ant64_bm3_ss25ka", "H", None, None, None, None, "inputs.json", 1, None, 1, "ant"),
    ("sico_specs_repo_ant64_b2_future09_ctrl.h", "repo_ant64_b2_future09_ctrl", "H", "5.e-2", None, None, None, "inputs.json", 1, None, 1, "ant"),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "H", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "H", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "H", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZSC.h", "grl40_bm5_paleo17a_CT4_BH0_ZSC", "H", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZSC.h", "grl40_bm5_paleo17a_CT4_BH0_FZSC", "H", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZLC.h", "grl40_bm5_paleo17a_CT4_BH0_ZLC", "H", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZLC.h", "grl40_bm5_paleo17a_CT4_BH0_FZLC", "H", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "H", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "H", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "H", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_BM5.h", "repo_grl16_bm5_ss25ka_BM5", "H", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZSC.h", "repo_grl16_bm5_ss25ka_ZSC", "H", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZLC.h", "repo_grl16_bm5_ss25ka_ZLC", "H", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "H", None, "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "H", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "H", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),

    # zs cases
    ("sico_specs_repo_ant64_bm3_ss25ka.h", "repo_ant64_bm3_ss25ka", "zs", None, None, None, None, "inputs.json", 1, None, 1, "ant"),
    ("sico_specs_repo_ant64_b2_future09_ctrl.h", "repo_ant64_b2_future09_ctrl", "zs", "5.e-2", None, None, None, "inputs.json", 1, None, 1, "ant"),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "zs", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "zs", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "zs", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZSC.h", "grl40_bm5_paleo17a_CT4_BH0_ZSC", "zs", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZSC.h", "grl40_bm5_paleo17a_CT4_BH0_FZSC", "zs", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZLC.h", "grl40_bm5_paleo17a_CT4_BH0_ZLC", "zs", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZLC.h", "grl40_bm5_paleo17a_CT4_BH0_FZLC", "zs", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "zs", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "zs", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "zs", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_BM5.h", "repo_grl16_bm5_ss25ka_BM5", "zs", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZSC.h", "repo_grl16_bm5_ss25ka_ZSC", "zs", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZLC.h", "repo_grl16_bm5_ss25ka_ZLC", "zs", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "zs", None, "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "zs", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "zs", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),

    # zl cases
    ("sico_specs_repo_ant64_bm3_ss25ka.h", "repo_ant64_bm3_ss25ka", "zl", None, None, None, None, "inputs.json", 1, None, 1, "ant"),
    ("sico_specs_repo_ant64_b2_future09_ctrl.h", "repo_ant64_b2_future09_ctrl", "zl", "5.e-2", None, None, None, "inputs.json", 1, None, 1, "ant"),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "zl", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "zl", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "zl", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZSC.h", "grl40_bm5_paleo17a_CT4_BH0_ZSC", "zl", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZSC.h", "grl40_bm5_paleo17a_CT4_BH0_FZSC", "zl", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZLC.h", "grl40_bm5_paleo17a_CT4_BH0_ZLC", "zl", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZLC.h", "grl40_bm5_paleo17a_CT4_BH0_FZLC", "zl", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "zl", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "zl", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "zl", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_BM5.h", "repo_grl16_bm5_ss25ka_BM5", "zl", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZSC.h", "repo_grl16_bm5_ss25ka_ZSC", "zl", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZLC.h", "repo_grl16_bm5_ss25ka_ZLC", "zl", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "zl", None, "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "zl", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "zl", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),

    # zl0 cases
    ("sico_specs_repo_ant64_bm3_ss25ka.h", "repo_ant64_bm3_ss25ka", "zl0", None, None, None, None, "inputs.json", 1, None, 1, "ant"),
    ("sico_specs_repo_ant64_b2_future09_ctrl.h", "repo_ant64_b2_future09_ctrl", "zl0", "5.e-2", None, None, None, "inputs.json", 1, None, 1, "ant"),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "zl0", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "zl0", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "zl0", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZSC.h", "grl40_bm5_paleo17a_CT4_BH0_ZSC", "zl0", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZSC.h", "grl40_bm5_paleo17a_CT4_BH0_FZSC", "zl0", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZLC.h", "grl40_bm5_paleo17a_CT4_BH0_ZLC", "zl0", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZLC.h", "grl40_bm5_paleo17a_CT4_BH0_FZLC", "zl0", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "zl0", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "zl0", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "zl0", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_BM5.h", "repo_grl16_bm5_ss25ka_BM5", "zl0", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZSC.h", "repo_grl16_bm5_ss25ka_ZSC", "zl0", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZLC.h", "repo_grl16_bm5_ss25ka_ZLC", "zl0", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "zl0", None, "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "zl0", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "zl0", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),

    # zb cases
    ("sico_specs_repo_ant64_bm3_ss25ka.h", "repo_ant64_bm3_ss25ka", "zb", None, None, None, None, "inputs.json", 1, None, 1, "ant"),
    ("sico_specs_repo_ant64_b2_future09_ctrl.h", "repo_ant64_b2_future09_ctrl", "zb", "5.e-2", None, None, None, "inputs.json", 1, None, 1, "ant"),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0.h", "grl40_bm5_paleo17a_CT4_BH0", "zb", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_BM5.h", "grl40_bm5_paleo17a_CT4_BH0_BM5", "zb", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FBM5.h", "grl40_bm5_paleo17a_CT4_BH0_FBM5", "zb", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZSC.h", "grl40_bm5_paleo17a_CT4_BH0_ZSC", "zb", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZSC.h", "grl40_bm5_paleo17a_CT4_BH0_FZSC", "zb", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_ZLC.h", "grl40_bm5_paleo17a_CT4_BH0_ZLC", "zb", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FZLC.h", "grl40_bm5_paleo17a_CT4_BH0_FZLC", "zb", "1.e-4", None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "zb", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "zb", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "zb", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_BM5.h", "repo_grl16_bm5_ss25ka_BM5", "zb", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZSC.h", "repo_grl16_bm5_ss25ka_ZSC", "zb", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZLC.h", "repo_grl16_bm5_ss25ka_ZLC", "zb", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "zb", None, "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_SVC.h", "grl40_bm5_paleo17a_CT4_BH0_SVC", "zb", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FSVC.h", "grl40_bm5_paleo17a_CT4_BH0_FSVC", "zb", "1.e-4", "0.01", None, None, "inputs.json", 1, None, 1, None),

    # q_geo cases
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "q_geo", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_BM5.h", "repo_grl16_bm5_ss25ka_BM5", "q_geo", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZSC.h", "repo_grl16_bm5_ss25ka_ZSC", "q_geo", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_ZLC.h", "repo_grl16_bm5_ss25ka_ZLC", "q_geo", None, None, None, None, "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "q_geo", None, "0.01", None, None, "inputs.json", 1, None, 1, None),

    # temp_c and age_c cases
    ("sico_specs_repo_grl16_bm5_ss25ka.h", "repo_grl16_bm5_ss25ka", "temp_c", None, None, "3", "40", "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_AC.h", "grl40_bm5_paleo17a_CT4_BH0_AC", "age_c", None, None, "3", "40", "inputs.json", 1, None, 1, None),
    ("sico_specs_grl40_bm5_paleo17a_CT4_BH0_FAC.h", "grl40_bm5_paleo17a_CT4_BH0_FAC", "age_c", None, None, "3", "40", "inputs.json", 1, None, 1, None),
    ("sico_specs_repo_grl16_bm5_ss25ka_CM3_AC.h", "repo_grl16_bm5_ss25ka_CM3_AC", "age_c", None, None, "3", "40", "inputs.json", 1, None, 1, None),
])
def test_tapenade_config(header_filename, head, iv, delta, low, dim, z, jsf, asi, lbfs, bia, dom):
    subprocess.run(['cp', f'headers/{header_filename}', '../headers'])
    subprocess.run(['cp', f'headers/{header_filename}', '../src/sico_specs.h'])

    cmd = f'python3 tapenade_config.py -head {head} -iv {iv} -jsf {jsf}'

    if delta is not None:
        cmd += f' -delta {delta}'
    if low is not None:
        cmd += f' -low {low}'
    if dim is not None:
        cmd += f' -dim {dim}'
    if z is not None:
        cmd += f' -z {z}'
    if lbfs is not None:
        cmd += f' -lbfs {lbfs}'
    if bia is not None:
        cmd += f' -bia {bia}'
    if asi is not None:
        cmd += f' -asi {asi}'
    if dom is not None:
        cmd += f' -dom {dom}'
    subprocess.run(cmd, shell=True, check=True)
