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

    ! Print matrices
    call print_matrix('D1MO', D1MO(1:6,1:6), ouf, dfmt)     ; write(ouf,*)
    call print_matrix('H1MO', H1MO(1:6,1:6), ouf, dfmt)     ; write(ouf,*)
    call print_matrix('D2MO', D2MO(1,1,1:6,1:6), ouf, dfmt) ; write(ouf,*)
    call print_matrix('H2MO', H2MO(1,1,1:6,1:6), ouf, dfmt) ; write(ouf,*)

    ! Obtain the natural occupation numbers: the eigenvalues of D^(1)
    call diag_2D_mat(D1MO, NONMO)
    call print_vector('Natural occupation numbers:', NONMO, 'row', ouf, '(*(f6.2))')
    write(ouf,*) ; write(ouf,'(a,f6.2)') 'Tr[D1] =', trace(D1MO) ; write(ouf,*)

    ! Compute the non-ee part of the energy: tr[D^(1) H^(1)]
    call H1D1_energy(D1MO, H1MO, EH1D1)
    write(ouf,'(a, d10.2, a)') 'EH1D1 = ', EH1D1, 'a.u.' ; write(ouf,*)

    ! Compute the total energy
    ! call energy(EH1D1, Eee_D2, ED2)

    ! Print the HF energy
    ! write(ouf,'(a, d10.2, a)') 'Exact energy = ', ED2, 'a.u.'

    ! Deallocate arrays
    deallocate(D1MO)
    deallocate(H1MO)
    deallocate(D2MO)
    deallocate(H2MO)


    !
    stop
endprogram main
