! +--------------------------------------------+
! | Author: Jose Antonio Quinonero Gris        |
! | Creation date: Tuesday 09:11:31 16/01/2024 |
! +--------------------------------------------+
module declarations 
    ! 
    ! Module to declare variables
    !
    implicit none
    logical :: tstdbg, diag
    integer :: i, j 
    integer :: termwidth, ouf, pmatsize
    character(len=80) :: d1fn, d2fn, h1fn, h2fn, dfmt
    real(kind=8) :: fthres 
    real(kind=8) :: EH1D1, EH2D2, E_exact
    real(kind=8) :: EeeELS, EeeBBC1, EeeBBC2, EeeBBC3, EeeBBC3M
    real(kind=8) :: EELS, EBBC1, EBBC2, EBBC3, EBBC3M
    real(kind=8), dimension(:), allocatable :: NONMO, NONSO 
    real(kind=8), dimension(:,:), allocatable :: D1MO, D1SO, H1MO, H1SO
    real(kind=8), dimension(:,:,:,:), allocatable :: D2MO, D2SO, H2MO, H2SO

    contains

    subroutine parameters
        implicit none

        ! Print debug information
        tstdbg = .true.
        
        ! floating-point numbers threshold to consider it 0, i.e.
        ! a -> 0 if a < fthres
        fthres = 1.d-10

        ! Width of the terminal/file to be written on
        termwidth = 46

        ! Logic unit of the output file
        ouf = 6

        ! D format
        dfmt = '(*(d10.2))'

        ! Matrix size to print: print A(1:pmatsize,1:pmatsize)
        pmatsize = 4

        return
    end subroutine parameters 

    !
end module declarations 
