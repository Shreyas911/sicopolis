  subroutine template()

    use OAD_tape
    use OAD_rev
     use oad_sico_variables_m, only : pi
!$TEMPLATE_PRAGMA_DECLARATIONS


    type(modeType) :: our_orig_mode
    real(8)              :: tempval
!
!       **** Statements ****
!


    integer iaddr
    external iaddr

    if (our_rev_mode%plain) then
! original function
      retval%v = erfc(x%v)
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
      if (oad_dt_sz.lt.oad_dt_ptr) call oad_dt_grow()
      oad_dt(oad_dt_ptr) = -1*2*exp(-1*(x%v*x%v))/sqrt(pi)
      oad_dt_ptr = oad_dt_ptr+1
      retval%v = erfc(x%v)
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
      oad_dt_ptr = oad_dt_ptr-1
      tempval = oad_dt(oad_dt_ptr)
      x%d = x%d+tempval*retval%d
      retval%d = 0.0
      print *, "After  sorsAD integer_tape_pointer ", oad_it_ptr
      print *, "After  sorsAD double_tape_pointer  ", oad_dt_ptr
      our_rev_mode=our_orig_mode
    end if
  end subroutine template
