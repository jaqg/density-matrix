! +--------------------------------------------+
! | Author: Jose Antonio Quinonero Gris        |
! | Creation date: Tuesday 17:27:38 30/04/2024 |
! +--------------------------------------------+
program main
    !
    ! Program to compute the HF energy from the 1e and 2e-density matrices
    !
    ! Modules
    use declarations
    use strings
    use IO
    use density_matrices
    !
    implicit none
    !
    ! === START OF THE PROGRAM ===
    !
     
    ! Initialize the parameters
    call parameters

    ! Print the title
    call title(ouf, 'Density Matrices', termwidth)

    ! Read the input
    call read_input('input', d1fn, d2fn, h1fn, h2fn)

    ! Read the density & integral matrices
    call read_2D_matrix(d1fn, .true., D1MO)  ! one-electron density matrix
    call read_2D_matrix(h1fn, .true., H1MO)  ! one-electron integral matrix
    call read_4D_matrix(d2fn, .false., D2MO)  ! two-electron density matrix
    call read_4D_matrix(h2fn, .false., H2MO)  ! test
    ! call read_4D_matrix(h2fn, .true., H2MO)  ! two-electron density matrix

    call print_matrix(D1MO(1:4,1:4), '(*(d10.2))') ; write(*,*)
    call print_matrix(H1MO(1:4,1:4), '(*(d10.2))') ; write(*,*)
    call print_matrix(D2MO(1,1,1:4,1:4), '(*(d10.2))') ; write(*,*)
    call print_matrix(H2MO(1,1,1:4,1:4), '(*(d10.2))') ; write(*,*)

    ! Deallocate arrays
    deallocate(D1MO)
    deallocate(H1MO)
    deallocate(D2MO)
    deallocate(H2MO)


    !
    stop
endprogram main
