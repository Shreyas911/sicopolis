  subroutine template()

    use OAD_tape
    use OAD_rev
    use sico_maths_m, only : tri_sle_local
    use sico_maths_m_grad, only : tri_sle_grad
!$TEMPLATE_PRAGMA_DECLARATIONS

    type(modeType) :: our_orig_mode
    integer iaddr
    external iaddr
    integer :: i, nrows_local
    if (our_rev_mode%plain) then
      call tri_sle_local(a0%v, a1%v, a2%v, x%v, b%v, nrows)
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
!!call push_i_s0(nrows)
      !integer_tape(integer_tape_pointer)=nrows
      !integer_tape_pointer=integer_tape_pointer+1
      if (oad_it_sz.lt.oad_it_ptr) call oad_it_grow()
      oad_it(oad_it_ptr) = nrows
      oad_it_ptr = oad_it_ptr+1
      do i=1,nrows
!call push_s0(a0(i)%v)
        !double_tape(double_tape_pointer)=a0(i)%v
        !double_tape_pointer=double_tape_pointer+1
        if (oad_dt_sz.lt.oad_dt_ptr) call oad_dt_grow()
        oad_dt(oad_dt_ptr) = a0(i)%v
        oad_dt_ptr = oad_dt_ptr+1
      end do
      do i=0,nrows
!call push_s0(a2(i)%v)
        !double_tape(double_tape_pointer)=a2(i)%v
        !double_tape_pointer=double_tape_pointer+1
        if (oad_dt_sz.lt.oad_dt_ptr) call oad_dt_grow()
        oad_dt(oad_dt_ptr) = a2(i)%v
        oad_dt_ptr = oad_dt_ptr+1
      end do
      call tri_sle_local(a0%v, a1%v, a2%v, x%v, b%v, nrows)
      do i=0,nrows
!call push_s0(a1(i)%v)
        !double_tape(double_tape_pointer)=a1(i)%v
        !double_tape_pointer=double_tape_pointer+1
        if (oad_dt_sz.lt.oad_dt_ptr) call oad_dt_grow()
        oad_dt(oad_dt_ptr) = a1(i)%v
        oad_dt_ptr = oad_dt_ptr+1
      end do
      do i=0,nrows
!call push_s0(b(i)%v)
        !double_tape(double_tape_pointer)=b(i)%v
        !double_tape_pointer=double_tape_pointer+1
        if (oad_dt_sz.lt.oad_dt_ptr) call oad_dt_grow()
        oad_dt(oad_dt_ptr) = b(i)%v
        oad_dt_ptr = oad_dt_ptr+1
      end do
      do i=0,nrows
!call push_s0(x(i)%v)
        !double_tape(double_tape_pointer)=x(i)%v
        !double_tape_pointer=double_tape_pointer+1
        if (oad_dt_sz.lt.oad_dt_ptr) call oad_dt_grow()
        oad_dt(oad_dt_ptr) = x(i)%v
        oad_dt_ptr = oad_dt_ptr+1
      end do
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
      !nrows_local=integer_tape(integer_tape_pointer)
      oad_it_ptr = oad_it_ptr-1
      nrows_local = oad_it(oad_it_ptr)
      do i=nrows_local, 0, -1
!call pop_s0(x(i)%v)
        !double_tape_pointer=double_tape_pointer-1
        !x(i)%v=double_tape(double_tape_pointer)
        oad_dt_ptr = oad_dt_ptr-1
        x(i)%v = oad_dt(oad_dt_ptr)
      end do
      do i=nrows_local, 0, -1
!call pop_s0(b(i)%v)
        !double_tape_pointer=double_tape_pointer-1
        !b(i)%v=double_tape(double_tape_pointer)
        oad_dt_ptr = oad_dt_ptr-1
        b(i)%v = oad_dt(oad_dt_ptr)
      end do
      do i=nrows_local, 0, -1
!call pop_s0(a1(i)%v)
        !double_tape_pointer=double_tape_pointer-1
        !a1(i)%v=double_tape(double_tape_pointer)
        oad_dt_ptr = oad_dt_ptr-1
        a1(i)%v = oad_dt(oad_dt_ptr)
      end do

      do i=nrows_local, 0, -1
!call pop_s0(a2(i)%v)
        !double_tape_pointer=double_tape_pointer-1
        !a2(i)%v=double_tape(double_tape_pointer)
        oad_dt_ptr = oad_dt_ptr-1
        a2(i)%v = oad_dt(oad_dt_ptr)
      end do
      do i=nrows_local, 1, -1
!call pop_s0(a0(i)%v)
        !double_tape_pointer=double_tape_pointer-1
        !a0(i)%v=double_tape(double_tape_pointer)
        oad_dt_ptr = oad_dt_ptr-1
        a0(i)%v = oad_dt(oad_dt_ptr)
      end do
      call tri_sle_grad(a0%v, a0%d, a1%v, a1%d, a2%v, a2%d, x%v, x%d, b%v, b%d, &
     &nrows_local)
      our_rev_mode=our_orig_mode
    end if
  end subroutine 
