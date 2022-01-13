module globals
        use sico_types_m
        implicit none
        real(dp), dimension(0:JMAX,0:IMAX) :: dzs_dx_aux, dzs_dy_aux

  integer(i4b) :: disc_DW
  integer(i4b) :: n_discharge_call_DW, iter_mar_coa_DW
  real(dp)     :: c_dis_0_DW, s_dis_DW, c_dis_fac_DW
  real(dp)     :: T_sub_PD_DW, alpha_sub_DW, alpha_o_DW, m_H_DW, m_D_DW, r_mar_eff_DW
  real(dp)     :: T_sea_freeze_DW
  real(dp)     :: dT_glann, dT_sub
  integer(i4b), dimension(0:JMAX,0:IMAX) :: mask_mar
  real(dp),     dimension(0:JMAX,0:IMAX) :: c_dis_DW
  real(dp),     dimension(0:JMAX,0:IMAX) :: cst_dist, cos_grad_tc
  real(dp),     dimension(0:JMAX,0:IMAX) :: dis_perp 

end module globals
