  subroutine template()

    use OAD_tape
    use OAD_rev
    use sico_maths_m, only : sor_sprs_local
    use sico_maths_m_grad, only : sor_sprs_grad
!$TEMPLATE_PRAGMA_DECLARATIONS


    type(modeType) :: our_orig_mode
INTEGER :: iter

INTEGER(4) :: nr, k, i, nmax_local, nnz_local, n_sprs_local
real(8)    :: omega_local, eps_sor_local
logical      :: flag_convergence
!
!       **** Statements ****
!


    integer iaddr
    external iaddr

    if (our_rev_mode%plain) then
! original function
      call sor_sprs_local(lgs_a_value%v, lgs_a_index, &
                    lgs_a_diag_index, lgs_a_ptr, &
                    lgs_b_value%v, nnz, nmax, n_sprs, omega, &
                    eps_sor, lgs_x_value%v, ierr)
    end if
    if (our_rev_mode%tape) then
! taping
! set up for plain execution
      our_orig_mode=our_rev_mode
      our_rev_mode%arg_store=.FALSE.
      our_rev_mode%arg_restore=.FALSE.
      our_rev_mode%plain=.TRUE.
      our_rev_mode%tape=.FALSE.
      our_rev_mode%adjoint=.FALSE.
      print *, "Before sorsTP integer_tape_pointer ", oad_it_ptr 
      print *, "Before sorsTP double_tape_pointer  ", oad_dt_ptr
      do i=1,nmax
        !double_tape(double_tape_pointer)=lgs_x_value(i)%v
        !double_tape_pointer=double_tape_pointer+1
        if (oad_dt_sz.lt.oad_dt_ptr) call oad_dt_grow()
        oad_dt(oad_dt_ptr) = lgs_x_value(i)%v
        oad_dt_ptr = oad_dt_ptr+1
      end do
      do i=1,nnz
        !double_tape(double_tape_pointer)=lgs_a_value(i)%v
        !double_tape_pointer=double_tape_pointer+1
        if (oad_dt_sz.lt.oad_dt_ptr) call oad_dt_grow()
        oad_dt(oad_dt_ptr) = lgs_a_value(i)%v
        oad_dt_ptr = oad_dt_ptr+1
      end do
      do i=1,nmax+1
        !integer_tape(integer_tape_pointer)=lgs_a_ptr(i)
        !integer_tape_pointer=integer_tape_pointer+1
        if (oad_it_sz.lt.oad_it_ptr) call oad_it_grow()
        oad_it(oad_it_ptr) = lgs_a_ptr(i)
        oad_it_ptr = oad_it_ptr+1
      end do
      do i=1,nnz
        !integer_tape(integer_tape_pointer)=lgs_a_index(i)
        !integer_tape_pointer=integer_tape_pointer+1
        if (oad_it_sz.lt.oad_it_ptr) call oad_it_grow()
        oad_it(oad_it_ptr) = lgs_a_index(i)
        oad_it_ptr = oad_it_ptr+1
      end do
      do i=1,nmax
        !integer_tape(integer_tape_pointer)=lgs_a_diag_index(i)
        !integer_tape_pointer=integer_tape_pointer+1
        if (oad_it_sz.lt.oad_it_ptr) call oad_it_grow()
        oad_it(oad_it_ptr) = lgs_a_diag_index(i)
        oad_it_ptr = oad_it_ptr+1
      end do
      do i=1,nmax
        !double_tape(double_tape_pointer)=lgs_b_value(i)%v
        !double_tape_pointer=double_tape_pointer+1
        if (oad_dt_sz.lt.oad_dt_ptr) call oad_dt_grow()
        oad_dt(oad_dt_ptr) = lgs_b_value(i)%v
        oad_dt_ptr = oad_dt_ptr+1
      end do
      call sor_sprs_local(lgs_a_value%v, lgs_a_index, &
                    lgs_a_diag_index, lgs_a_ptr, &
                    lgs_b_value%v, nnz, nmax, n_sprs, omega, &
                    eps_sor, lgs_x_value%v, ierr)
      do i=1,nmax
        !double_tape(double_tape_pointer)=lgs_x_value(i)%v
        !double_tape_pointer=double_tape_pointer+1
        if (oad_dt_sz.lt.oad_dt_ptr) call oad_dt_grow()
        oad_dt(oad_dt_ptr) = lgs_x_value(i)%v
        oad_dt_ptr = oad_dt_ptr+1
      end do
      !double_tape(double_tape_pointer)=omega
      !double_tape_pointer=double_tape_pointer+1
        if (oad_dt_sz.lt.oad_dt_ptr) call oad_dt_grow()
        oad_dt(oad_dt_ptr) = omega
        oad_dt_ptr = oad_dt_ptr+1
      !double_tape(double_tape_pointer)=eps_sor
      !double_tape_pointer=double_tape_pointer+1
        if (oad_dt_sz.lt.oad_dt_ptr) call oad_dt_grow()
        oad_dt(oad_dt_ptr) = eps_sor
        oad_dt_ptr = oad_dt_ptr+1
      !integer_tape(integer_tape_pointer)=n_sprs
      !integer_tape_pointer=integer_tape_pointer+1
        if (oad_it_sz.lt.oad_it_ptr) call oad_it_grow()
        oad_it(oad_it_ptr) = n_sprs
        oad_it_ptr = oad_it_ptr+1
      !integer_tape(integer_tape_pointer)=nnz
      !integer_tape_pointer=integer_tape_pointer+1
        if (oad_it_sz.lt.oad_it_ptr) call oad_it_grow()
        oad_it(oad_it_ptr) = nnz
        oad_it_ptr = oad_it_ptr+1
      !integer_tape(integer_tape_pointer)=nmax
      !integer_tape_pointer=integer_tape_pointer+1
        if (oad_it_sz.lt.oad_it_ptr) call oad_it_grow()
        oad_it(oad_it_ptr) = nmax
        oad_it_ptr = oad_it_ptr+1
      our_rev_mode=our_orig_mode
    end if
    if (our_rev_mode%adjoint) then
! adjoint
      our_orig_mode=our_rev_mode
      our_rev_mode%arg_store=.FALSE.
      our_rev_mode%arg_restore=.FALSE.
      our_rev_mode%plain=.TRUE.
      our_rev_mode%tape=.FALSE.
      our_rev_mode%adjoint=.FALSE.
      !integer_tape_pointer=integer_tape_pointer-1
      !nmax_local=integer_tape(integer_tape_pointer)
      oad_it_ptr = oad_it_ptr-1
      nmax_local = oad_it(oad_it_ptr)
      !integer_tape_pointer=integer_tape_pointer-1
      !nnz_local=integer_tape(integer_tape_pointer)
      oad_it_ptr = oad_it_ptr-1
      nnz_local = oad_it(oad_it_ptr)
      !integer_tape_pointer=integer_tape_pointer-1
      !n_sprs_local=integer_tape(integer_tape_pointer)
      oad_it_ptr = oad_it_ptr-1
      n_sprs_local = oad_it(oad_it_ptr)
      !double_tape_pointer=double_tape_pointer-1
      !eps_sor_local=double_tape(double_tape_pointer)
      oad_dt_ptr = oad_dt_ptr-1
      eps_sor_local = oad_dt(oad_dt_ptr)
      !double_tape_pointer=double_tape_pointer-1
      !omega_local=double_tape(double_tape_pointer)
      oad_dt_ptr = oad_dt_ptr-1
      omega_local = oad_dt(oad_dt_ptr)
      do i=nmax_local, 1, -1
        !double_tape_pointer=double_tape_pointer-1
        !lgs_x_value(i)%v=double_tape(double_tape_pointer)
        oad_dt_ptr = oad_dt_ptr-1
        lgs_x_value(i)%v = oad_dt(oad_dt_ptr)
      end do
      do i=nmax_local, 1, -1
        !double_tape_pointer=double_tape_pointer-1
        !lgs_b_value(i)%v=double_tape(double_tape_pointer)
        oad_dt_ptr = oad_dt_ptr-1
        lgs_b_value(i)%v = oad_dt(oad_dt_ptr)
      end do
      do i=nmax_local, 1, -1
        !integer_tape_pointer=integer_tape_pointer-1
        !lgs_a_diag_index(i)=integer_tape(integer_tape_pointer)
        oad_it_ptr = oad_it_ptr-1
        lgs_a_diag_index(i) = oad_it(oad_it_ptr)
      end do
      do i=nnz_local, 1, -1
        !integer_tape_pointer=integer_tape_pointer-1
        !lgs_a_index(i)=integer_tape(integer_tape_pointer)
        oad_it_ptr = oad_it_ptr-1
        lgs_a_index(i) = oad_it(oad_it_ptr)
      end do
      do i=nmax_local+1, 1, -1
        !integer_tape_pointer=integer_tape_pointer-1
        !lgs_a_ptr(i)=integer_tape(integer_tape_pointer)
        oad_it_ptr = oad_it_ptr-1
        lgs_a_ptr(i) = oad_it(oad_it_ptr)
      end do
      do i=nnz_local, 1, -1
        !double_tape_pointer=double_tape_pointer-1
        !lgs_a_value(i)%v=double_tape(double_tape_pointer)
        oad_dt_ptr = oad_dt_ptr-1
        lgs_a_value(i)%v = oad_dt(oad_dt_ptr)
      end do
      call sor_sprs_grad(lgs_a_value%v, lgs_a_value%d,&
          lgs_a_index, lgs_a_diag_index, lgs_a_ptr, &
          lgs_b_value%v, lgs_b_value%d, &
          nnz_local, nmax_local, n_sprs_local,&
           omega_local, eps_sor_local, &
          lgs_x_value%v, lgs_x_value%d, ierr)
      do i=nmax_local, 1, -1
        !double_tape_pointer=double_tape_pointer-1
        !lgs_x_value(i)%v=double_tape(double_tape_pointer)
        oad_dt_ptr = oad_dt_ptr-1
        lgs_x_value(i)%v = oad_dt(oad_dt_ptr)
      end do
      print *, "After  sorsAD integer_tape_pointer ", oad_it_ptr
      print *, "After  sorsAD double_tape_pointer  ", oad_dt_ptr
      our_rev_mode=our_orig_mode
    end if
  end subroutine template
