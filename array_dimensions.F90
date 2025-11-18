module array_dimensions
    implicit none

    ! This module defines array dimensions at compile time
    ! kLength is set via preprocessor directive -DkLength=<value>
    ! Example: gfortran -DkLength=144 ...

    integer, parameter :: dp = 8
    integer, parameter :: KLENGTH_PARAM = KLENGTH_SIZE

end module array_dimensions
