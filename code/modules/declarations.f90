! +--------------------------------------------+
! | Author: Jose Antonio Quinonero Gris        |
! | Creation date: Tuesday 09:11:31 16/01/2024 |
! +--------------------------------------------+
module declarations 
    ! 
    ! Module to declare variables
    !
    implicit none
    integer :: i, j 
    integer :: termwidth, ouf
    character(len=80) :: d1fn, d2fn, h1fn, h2fn, dfmt
    real(kind=8) :: fthres 
    real(kind=8) :: EH1D1, EH2D2, E_exact
    real(kind=8) :: EeeBBC1, EeeBBC3, EBBC1, EBBC3 
    real(kind=8), dimension(:), allocatable :: NONMO, NONSO 
    real(kind=8), dimension(:,:), allocatable :: D1MO, H1MO, D1SO
    real(kind=8), dimension(:,:,:,:), allocatable :: D2MO, H2MO, H2SO

    contains

    subroutine parameters
        implicit none
        
        ! floating-point numbers threshold to consider it 0, i.e.
        ! a -> 0 if a < fthres
        fthres = 1.d-13

        ! Width of the terminal/file to be written on
        termwidth = 46

        ! Logic unit of the output file
        ouf = 6

        ! D format
        dfmt = '(*(d10.2))'

        return
    end subroutine parameters 

    !
end module declarations 
