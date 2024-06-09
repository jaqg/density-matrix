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

    ! Obtain the MO's occupation numbers: the eigenvalues of D^(1)
    call diag_2D_mat(D1MO, NONMO)
    call print_vector('MO occupation numbers:', NONMO, 'row', ouf, '(*(f6.2))')
    write(ouf,*)
    write(ouf,'(a,f6.2)') 'Tr[D1MO] =', trace(D1MO) ; write(ouf,*)

    ! Transform D^(1) from MO to spin-orbital basis
    call MO_to_SO(D1MO, D1SO)
    call print_matrix('D1SO', D1SO(1:6,1:6), ouf, dfmt)     ; write(ouf,*)

    ! Obtain the natural occupation numbers: the eigenvalues of D^(1) in SO
    call diag_2D_mat(D1SO, NONSO)
    call print_vector('natural occupation numbers:',NONSO,'row',ouf,'(*(f6.2))')
    write(ouf,*)
    write(ouf,'(a,f6.2)') 'Tr[D1SO] =', trace(D1SO) ; write(ouf,*)

    ! Compute the non-ee part of the energy: tr[D^(1) H^(1)]
    call H1D1_energy(D1MO, H1MO, EH1D1)
!   write(ouf,'(a, d10.2, a)') 'EH1D1 = ', EH1D1, 'a.u.' ; write(ouf,*)

    !
    ! Compute the Eee part with:
    !
    !   - Exact from the D^(2): EeeD2 = 1/2 sum_{pqrs} (pq|rs) D2_{pqrs}
    call Eee_D2(D2MO, H2MO, EH2D2)
    !   - BB functionals: BBC1, BBC3
    ! call Eee_BBC('BBC1', NONSO, H2SO, EeeBBC1)  ! BUENO
    ! call Eee_BBC('BBC3', NONSO, H2SO, EeeBBC3)  ! BUENO
    ! TODO: pasar H2MO -> H2SO
    call Eee_BBC('BBC1', NONMO, H2MO, EeeBBC1)  ! TEST
    call Eee_BBC('BBC3', NONMO, H2MO, EeeBBC3)  ! TEST


    ! Compute the total energy for each functional
    call energy(EH1D1, EH2D2, E_exact)
    call energy(EH1D1, EeeBBC1, EBBC1)
    call energy(EH1D1, EeeBBC3, EBBC3)

    ! Print the results
    write(ouf,'(a, d10.2, a)') 'Exact energy = ', E_exact, 'a.u.'
    write(ouf,'(a, d10.2, a)') 'BBC1 energy = ', EBBC1, 'a.u.'
    write(ouf,'(a, d10.2, a)') 'BBC3 energy = ', EBBC3, 'a.u.'

    ! Deallocate arrays
    deallocate(D1MO)
    deallocate(H1MO)
    deallocate(D2MO)
    deallocate(H2MO)


    !
    stop
endprogram main
