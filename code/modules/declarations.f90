! +--------------------------------------------+
! | Author: Jose Antonio Quinonero Gris        |
! | Creation date: Tuesday 09:11:31 16/01/2024 |
! +--------------------------------------------+
module declarations 
    ! 
    ! Module to declare variables
    !
    implicit none
    integer :: termwidth, ouf
    character(len=80) :: d1fn, d2fn, h1fn, h2fn 
    real(kind=8), dimension(:,:), allocatable :: D1MO, H1MO
    real(kind=8), dimension(:,:,:,:), allocatable :: D2MO, H2MO

    contains

    subroutine parameters
        implicit none
        
        ! Width of the terminal/file to be written on
        termwidth = 46

        ! Logic unit of the output file
        ouf = 6

        return
    end subroutine parameters 

    !
end module declarations 
