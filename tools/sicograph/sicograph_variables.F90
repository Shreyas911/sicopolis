!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!   Module     :  s i c o g r a p h _ v a r i a b l e s

!   Purpose    :  Declarations of global variables for SICOGRAPH.

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module sicograph_variables

use sicograph_types

integer(i1b) :: maske_gr(0:IMAX,0:JMAX), n_cts_gr(0:IMAX,0:JMAX)
integer(i4b) :: kc_cts_gr(0:IMAX,0:JMAX)
real(sp) :: time_gr, delta_ts_gr, glac_index_gr, z_sl_gr, &
            V_tot_gr, A_grounded_gr, A_floating_gr, &
            xi_gr(0:IMAX), eta_gr(0:JMAX), &
            z_c_gr(0:IMAX,0:JMAX,0:KCMAX), &
            z_t_gr(0:IMAX,0:JMAX,0:KTMAX), &
            z_r_gr(0:IMAX,0:JMAX,0:KRMAX), &
            zs_gr(0:IMAX,0:JMAX), zm_gr(0:IMAX,0:JMAX), &
            zb_gr(0:IMAX,0:JMAX), zl_gr(0:IMAX,0:JMAX), &
            H_cold_gr(0:IMAX,0:JMAX), H_temp_gr(0:IMAX,0:JMAX), &
            H_gr(0:IMAX,0:JMAX), H_R_gr, &
            zl0_gr(0:IMAX,0:JMAX), &
            vx_c_gr(0:IMAX,0:JMAX,0:KCMAX), &
            vy_c_gr(0:IMAX,0:JMAX,0:KCMAX), &
            vz_c_gr(0:IMAX,0:JMAX,0:KCMAX), &
            vx_t_gr(0:IMAX,0:JMAX,0:KTMAX), &
            vy_t_gr(0:IMAX,0:JMAX,0:KTMAX), &
            vz_t_gr(0:IMAX,0:JMAX,0:KTMAX), &
            temp_c_gr(0:IMAX,0:JMAX,0:KCMAX), &
            temph_c_gr(0:IMAX,0:JMAX,0:KCMAX), &
            age_c_gr(0:IMAX,0:JMAX,0:KCMAX), &
            enth_c_gr(0:IMAX,0:JMAX,0:KCMAX), &
            omega_c_gr(0:IMAX,0:JMAX,0:KCMAX), &
            enh_c_gr(0:IMAX,0:JMAX,0:KCMAX), &
            enth_t_gr(0:IMAX,0:JMAX,0:KTMAX), &
            omega_t_gr(0:IMAX,0:JMAX,0:KTMAX), &
            temp_t_gr(0:IMAX,0:JMAX,0:KTMAX), &
            age_t_gr(0:IMAX,0:JMAX,0:KTMAX), &
            enh_t_gr(0:IMAX,0:JMAX,0:KTMAX), &
            temp_r_gr(0:IMAX,0:JMAX,0:KRMAX), &
            vx_b_g_gr(0:IMAX,0:JMAX), vy_b_g_gr(0:IMAX,0:JMAX), &
            vz_b_gr(0:IMAX,0:JMAX), vh_b_gr(0:IMAX,0:JMAX), &
            vx_s_g_gr(0:IMAX,0:JMAX), vy_s_g_gr(0:IMAX,0:JMAX), &
            vz_s_gr(0:IMAX,0:JMAX), vh_s_gr(0:IMAX,0:JMAX), &
            temp_b_gr(0:IMAX,0:JMAX), temph_b_gr(0:IMAX,0:JMAX), &
            qx_gr(0:IMAX,0:JMAX), qy_gr(0:IMAX,0:JMAX), &
            Q_bm_gr(0:IMAX,0:JMAX), Q_tld_gr(0:IMAX,0:JMAX), &
            am_perp_gr(0:IMAX,0:JMAX), &
            dzs_dtau_gr(0:IMAX,0:JMAX), dzm_dtau_gr(0:IMAX,0:JMAX), &
            dzb_dtau_gr(0:IMAX,0:JMAX), dzl_dtau_gr(0:IMAX,0:JMAX), &
            dH_c_dtau_gr(0:IMAX,0:JMAX), dH_t_dtau_gr(0:IMAX,0:JMAX), &
            dH_dtau_gr(0:IMAX,0:JMAX), &
            precip_gr(0:IMAX,0:JMAX), precip_lgm_anom_gr(0:IMAX,0:JMAX), &
            temp_s_gr(0:IMAX,0:JMAX), temp_s_lgm_anom_gr(0:IMAX,0:JMAX), &
            as_perp_gr(0:IMAX,0:JMAX), qgeo_gr(0:IMAX,0:JMAX), &
            p_b_w_gr(0:IMAX,0:JMAX), H_w_gr(0:IMAX,0:JMAX), &
            q_gl_g_gr(0:IMAX,0:JMAX)
character (len=256) :: ch_attr_title, ch_attr_institution, ch_attr_source, &
                       ch_attr_history, ch_attr_references
logical :: read_erg_nc_eof_flag

real(sp) :: time(1048576), D_Ts(1048576), z_sl(1048576), &
            V_tot(1048576), V_grounded(1048576), V_floating(1048576), &
            A_tot(1048576), A_grounded(1048576), A_floating(1048576), &
            V_sle(1048576), V_temp(1048576), A_temp(1048576), &
            H_max(1048576), H_t_max(1048576), zs_max(1048576), &
            vs_max(1048576), Tbh_max(1048576)

#if defined(HEINO)
real(sp) :: H_ave_sed(1048576), Tbh_ave_sed(1048576), Atb_sed(1048576)
#endif

#if defined(ANT)
real(sp) :: H_Vo(1048576), vs_Vo(1048576), Tb_Vo(1048576), Rb_Vo(1048576), &
            H_DA(1048576), vs_DA(1048576), Tb_DA(1048576), Rb_DA(1048576), &
            H_DC(1048576), vs_DC(1048576), Tb_DC(1048576), Rb_DC(1048576), &
            H_DF(1048576), vs_DF(1048576), Tb_DF(1048576), Rb_DF(1048576), &
            H_Ko(1048576), vs_Ko(1048576), Tb_Ko(1048576), Rb_Ko(1048576), &
            H_By(1048576), vs_By(1048576), Tb_By(1048576), Rb_By(1048576)
#elif defined(EMTP2SGE)
real(sp) :: H_P1(1048576), vs_P1(1048576), Tb_P1(1048576), Rb_P1(1048576), &
            H_P2(1048576), vs_P2(1048576), Tb_P2(1048576), Rb_P2(1048576)
#elif defined(GRL)
real(sp) :: HHGR(1048576), vs_GR(1048576), Tb_GR(1048576), Rb_GR(1048576), &
            H_G2(1048576), vs_G2(1048576), Tb_G2(1048576), Rb_G2(1048576), &
            H_D3(1048576), vs_D3(1048576), Tb_D3(1048576), Rb_D3(1048576), &
            H_CC(1048576), vs_CC(1048576), Tb_CC(1048576), Rb_CC(1048576), &
            H_NG(1048576), vs_NG(1048576), Tb_NG(1048576), Rb_NG(1048576)
#elif defined(HEINO)
real(sp) :: H_P1(1048576), vs_P1(1048576), Tb_P1(1048576), Rb_P1(1048576), &
            H_P2(1048576), vs_P2(1048576), Tb_P2(1048576), Rb_P2(1048576), &
            H_P3(1048576), vs_P3(1048576), Tb_P3(1048576), Rb_P3(1048576), &
            H_P4(1048576), vs_P4(1048576), Tb_P4(1048576), Rb_P4(1048576), &
            H_P5(1048576), vs_P5(1048576), Tb_P5(1048576), Rb_P5(1048576), &
            H_P6(1048576), vs_P6(1048576), Tb_P6(1048576), Rb_P6(1048576), &
            H_P7(1048576), vs_P7(1048576), Tb_P7(1048576), Rb_P7(1048576)
#elif defined(NMARS)
real(sp) :: H_NP(1048576), vs_NP(1048576), Tb_NP(1048576), Rb_NP(1048576), &
            H_C1(1048576), vs_C1(1048576), Tb_C1(1048576), Rb_C1(1048576), &
            H_C2(1048576), vs_C2(1048576), Tb_C2(1048576), Rb_C2(1048576)
#endif

real(sp) :: ea, zeta_c(0:KCMAX), eaz_c(0:KCMAX), &
            eaz_c_quotient(0:KCMAX), zeta_t(0:KTMAX), zeta_r(0:KRMAX)

end module sicograph_variables
!
