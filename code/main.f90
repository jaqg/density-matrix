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
    use matmod
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
    if (tstdbg) then 
        call print_matrix('D1MO', D1MO(1:6,1:6), ouf, dfmt)    
        call print_matrix('H1MO', H1MO(1:6,1:6), ouf, dfmt)    
        call print_matrix('D2MO', D2MO(1,1,1:6,1:6), ouf, dfmt)
        call print_matrix('H2MO', H2MO(1,1,1:6,1:6), ouf, dfmt)
    end if

    ! Obtain the MO's occupation numbers: the eigenvalues of D^(1)
    call is_diag(D1MO, fthres, diag)  ! Check if D1MO is already diagonal
    if (diag) then
        call extract_diag(D1MO, NONMO)
    else
        write(ouf,'(a)') 'D1MO is not diagonal. Diagonalizing...'
        call diag_2D_mat(D1MO, .true., NONMO)
    end if
    if (tstdbg) then 
        call print_vector('MO occupation numbers:', NONMO, 'row', ouf, '(*(f6.2))')
        write(ouf,'(a)') 'Occupied MO orbitals:'
        j = 1
        do i = 1, size(NONMO)
            if (NONMO(i).gt.fthres) then
                write(ouf,'(2(i0,2x),f4.2)') j, i, NONMO(i)
                j = j + 1
            end if
        end do
        write(ouf,*)
        write(ouf,'(a,f6.2)') 'Tr[D1MO] =', trace(D1MO) ; write(ouf,*)
    end if
    ! TODO: implementar que imprima Nelec desde DALTON para poder usarlo aqui
    call checknorm_D1(D1MO, 10, fthres, ouf)

    ! Transform D^(1), D^(2), H1, H2 from MO to spin-orbital basis
    ! TODO: comprobar que la transformacion para H1 y H2 este bien
    call MO_to_SO_D1(D1MO, D1SO)
    call MO_to_SO_D2(D2MO, D2SO)
    call MO_to_SO_D1(H1MO, H1SO)
    call MO_to_SO_D2(H2MO, H2SO)

    if (tstdbg) then 
        call print_matrix('D1SO', D1SO(1:6,1:6), ouf, dfmt) 
        call print_matrix('D2SO', D2SO(1,1,1:6,1:6), ouf, dfmt) 
        call print_matrix('H1SO', H1SO(1:6,1:6), ouf, dfmt) 
        call print_matrix('H2SO', H2SO(1,1,1:6,1:6), ouf, dfmt) 
    end if

    ! Obtain the natural occupation numbers: the eigenvalues of D^(1) in SO
    call is_diag(D1SO, fthres, diag)
    if (diag) then
        call extract_diag(D1SO, NONSO)
    else
        write(ouf,'(a)') 'D1SO is not diagonal. Diagonalizing...'
        call diag_2D_mat(D1SO, .true., NONSO)
    end if
    if (tstdbg) then 
        call print_vector('natural occupation numbers:',NONSO,'row',ouf,'(*(f6.2))')
        write(ouf,'(a)') 'Occupied SO orbitals:'
        j = 1
        do i = 1, size(NONSO)
            if (NONSO(i).gt.fthres) then
                write(ouf,'(2(i0,2x),f4.2)') j, i, NONSO(i)
                j = j + 1
            end if
        end do
        write(ouf,*)
        write(ouf,'(a,f6.2)') 'Tr[D1SO] =', trace(D1SO) ; write(ouf,*)
    end if
    ! TODO: implementar que imprima Nelec desde DALTON para poder usarlo aqui
    call checknorm_D1(D1SO, 10, fthres, ouf)

    ! Compute the non-ee part of the energy: tr[D^(1) H^(1)]
    call H1D1_energy(D1MO, H1MO, EH1D1)
    if (tstdbg) then 
          write(ouf,'(a, d10.2, a)') 'EH1D1 = ', EH1D1, 'a.u.' ; write(ouf,*)
    end if

    !
    ! Compute the Eee part with:
    !
    !   - Exact from the D^(2): EeeD2 = 1/2 sum_{pqrs} (pq|rs) D2_{pqrs}
    call Eee_D2(D2MO, H2MO, EH2D2)
    !   - Extended LS functional: ELS
    ! TODO: pasar H2MO -> H2SO
    ! TODO: comprobar que H2SO y D2SO estan bien
    call Eee_ELS(NONSO, H2SO, EeeELS)  ! BUENO
    ! call Eee_ELS(NONMO, H2MO, EeeELS)  ! TEST
    !   - BB functionals: BBC1, BBC3
    ! call Eee_BBC('BBC1', NONSO, H2SO, EeeBBC1)  ! BUENO
    ! call Eee_BBC('BBC2', NONSO, H2SO, EeeBBC2)  ! BUENO
    ! call Eee_BBC('BBC3', NONSO, H2SO, EeeBBC3)  ! BUENO
    ! call Eee_BBC('BBC3M', NONSO, H2SO, EeeBBC3M)  ! BUENO
    call Eee_BBC('BBC1', NONMO/2.d0, H2MO, EeeBBC1)  ! TEST
    call Eee_BBC('BBC2', NONMO/2.d0, H2MO, EeeBBC2)  ! TEST
    call Eee_BBC('BBC3', NONMO/2.d0, H2MO, EeeBBC3)  ! TEST
    call Eee_BBC('BBC3M', NONMO/2.d0, H2MO, EeeBBC3M)  ! TEST
    ! call Eee_BBC('BBC1', NONMO, H2MO, EeeBBC1)  ! TEST
    ! call Eee_BBC('BBC2', NONMO, H2MO, EeeBBC2)  ! TEST
    ! call Eee_BBC('BBC3', NONMO, H2MO, EeeBBC3)  ! TEST
    ! call Eee_BBC('BBC3M', NONMO, H2MO, EeeBBC3M)  ! TEST
    !   - PNOF functionals: PNOF5
    ! call Eee_PNOF(NONSO, 1.d-16)


    ! Compute the total energy for each functional
    call energy(EH1D1, EH2D2, E_exact)
    call energy(EH1D1, EeeELS, EELS)
    call energy(EH1D1, EeeBBC1, EBBC1)
    call energy(EH1D1, EeeBBC2, EBBC2)
    call energy(EH1D1, EeeBBC3, EBBC3)
    call energy(EH1D1, EeeBBC3M, EBBC3M)

    ! Print the results
    if (tstdbg) then 
        write(ouf,'(a, f10.4, 1x, a)') 'Exact energy = ', E_exact, 'a.u.'
        write(ouf,'(a, f10.4, 1x, a)') 'ELS energy = ', EELS, 'a.u.'
        write(ouf,'(a, f10.4, 1x, a)') 'BBC1 energy = ', EBBC1, 'a.u.'
        write(ouf,'(a, f10.4, 1x, a)') 'BBC2 energy = ', EBBC2, 'a.u.'
        write(ouf,'(a, f10.4, 1x, a)') 'BBC3 energy = ', EBBC3, 'a.u.'
        write(ouf,'(a, f10.4, 1x, a)') 'BBC3M energy = ', EBBC3M, 'a.u.'
    end if

    ! Compute the Minkowski distances
    ! notation: MD2_D2SO: Minkowski Distance (order) 2 for D_SO^(2)
    call Minkowski_distance_4D(2, D2SO, D2SO, MD2_D2SO)


    ! Print the results
    call print_energies

    ! Deallocate arrays
    deallocate(D1MO)
    deallocate(H1MO)
    deallocate(D2MO)
    deallocate(H2MO)
    !
    stop
endprogram main
