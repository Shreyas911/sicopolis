!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!   Module     :  s i c o g r a p h _ t y p e s

!   Purpose    :  Declarations of kind types for SICOGRAPH.

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module sicograph_types

integer, parameter :: i1b = selected_int_kind(2)   ! 1-byte integers
!!! integer, parameter :: i2b = selected_int_kind(4)   ! 2-byte integers
integer, parameter :: i4b = selected_int_kind(9)   ! 4-byte integers
integer, parameter :: sp  = kind(1.0)              ! single-precision reals
integer, parameter :: dp  = kind(1.0d0)            ! double-precision reals

end module sicograph_types
!
