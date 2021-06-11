  subroutine template()

    use OAD_tape
    use OAD_rev
    use sico_maths_m, only : sico_lis_solver_local
    use sico_maths_m_grad, only : sico_lis_solver_grad
!$TEMPLATE_PRAGMA_DECLARATIONS

    type(modeType) :: our_orig_mode
INTEGER(4) :: nr, k, i, nmax_local, nnz_local
!
!       **** Statements ****
!


    integer iaddr
    external iaddr

    if (our_rev_mode%plain) then
! original function
      call sico_lis_solver_local(nmax, nnz, &
                           lgs_a_ptr, lgs_a_index, &
                           lgs_a_value%v, lgs_b_value%v, lgs_x_value%v)
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
      print *, "Before sico_lis_solverTP integer_tape_pointer ", oad_it_ptr
      print *, "Before sico_lis_solverTP double_tape_pointer  ", oad_dt_ptr
      do i=1,nmax
        if (oad_dt_sz.lt.oad_dt_ptr) call oad_dt_grow()
        oad_dt(oad_dt_ptr) = lgs_x_value(i)%v
        oad_dt_ptr = oad_dt_ptr+1
      end do
      do i=1,nnz
        if (oad_dt_sz.lt.oad_dt_ptr) call oad_dt_grow()
        oad_dt(oad_dt_ptr) = lgs_a_value(i)%v
        oad_dt_ptr = oad_dt_ptr+1
      end do
      do i=1,nmax+1
        if (oad_it_sz.lt.oad_it_ptr) call oad_it_grow()
        oad_it(oad_it_ptr) = lgs_a_ptr(i)
        oad_it_ptr = oad_it_ptr+1
      end do
      do i=1,nnz
        if (oad_it_sz.lt.oad_it_ptr) call oad_it_grow()
        oad_it(oad_it_ptr) = lgs_a_index(i)
        oad_it_ptr = oad_it_ptr+1
      end do
      do i=1,nmax
        if (oad_dt_sz.lt.oad_dt_ptr) call oad_dt_grow()
        oad_dt(oad_dt_ptr) = lgs_b_value(i)%v
        oad_dt_ptr = oad_dt_ptr+1
      end do
      call sico_lis_solver_local(nmax, nnz, &
                           lgs_a_ptr, lgs_a_index, &
                           lgs_a_value%v, lgs_b_value%v, lgs_x_value%v)
      !do i=1,nmax
      !  if (oad_dt_sz.lt.oad_dt_ptr) call oad_dt_grow()
      !  oad_dt(oad_dt_ptr) = lgs_x_value(i)%v
      !  oad_dt_ptr = oad_dt_ptr+1
      !end do
        if (oad_it_sz.lt.oad_it_ptr) call oad_it_grow()
        oad_it(oad_it_ptr) = nnz
        oad_it_ptr = oad_it_ptr+1
        if (oad_it_sz.lt.oad_it_ptr) call oad_it_grow()
        oad_it(oad_it_ptr) = nmax
        oad_it_ptr = oad_it_ptr+1
        write (*,*) "nmax", nmax, nnz
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

      oad_it_ptr = oad_it_ptr-1
      nmax_local = oad_it(oad_it_ptr)
      oad_it_ptr = oad_it_ptr-1
      nnz_local = oad_it(oad_it_ptr)
       write (*,*) "nmax", nmax_local, nnz_local
      !do i=nmax_local, 1, -1
      !  oad_dt_ptr = oad_dt_ptr-1
      !  lgs_x_value(i)%v = oad_dt(oad_dt_ptr)
      !end do
      do i=nmax_local, 1, -1
        oad_dt_ptr = oad_dt_ptr-1
        lgs_b_value(i)%v = oad_dt(oad_dt_ptr)
      end do
      do i=nnz_local, 1, -1
        oad_it_ptr = oad_it_ptr-1
        lgs_a_index(i) = oad_it(oad_it_ptr)
      end do
      do i=nmax_local+1, 1, -1
        oad_it_ptr = oad_it_ptr-1
        lgs_a_ptr(i) = oad_it(oad_it_ptr)
      end do
      do i=nnz_local, 1, -1
        oad_dt_ptr = oad_dt_ptr-1
        lgs_a_value(i)%v = oad_dt(oad_dt_ptr)
      end do
      do i=nmax_local, 1, -1
        oad_dt_ptr = oad_dt_ptr-1
        lgs_x_value(i)%v = oad_dt(oad_dt_ptr)
      end do
      call sico_lis_solver_grad(nmax_local, nnz_local, &
                           lgs_a_ptr, lgs_a_index, &
                           lgs_a_value%v, lgs_a_value%d, lgs_b_value%v,&
                           lgs_b_value%d, lgs_x_value%v, lgs_x_value%d)
      !do i=nmax_local, 1, -1
      !  oad_dt_ptr = oad_dt_ptr-1
      !  lgs_x_value(i)%v = oad_dt(oad_dt_ptr)
      !end do
      print *, "After sico_lis_solverTP integer_tape_pointer ", oad_it_ptr
      print *, "After sico_lis_solverTP double_tape_pointer  ", oad_dt_ptr
      our_rev_mode=our_orig_mode
    end if
  end subroutine template
