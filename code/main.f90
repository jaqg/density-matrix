program RDMFT
    
   ! Modules
   use declarations
   use strings
   use IO
   use matmod
   use density_matrices
   use LS
   use BBC

   implicit none

   ! +-----------------------------------------------------------------------+
   ! |                         START OF THE PROGRAM                          |
   ! +-----------------------------------------------------------------------+

   ! Initialize the parameters
   call parameters

   ! Print the title
   call title(ouf, 'Density Matrices', termwidth)

   ! Read the input
   write(ouf,'(a)') 'Reading input...'; write(ouf,*)
   call read_input('INPUT', sirifn, intfn)

   ! Open the energies output file
   call qopen('ENERGIES', 'formatted', 'sequential', eolu)
   write(eolu,'(11x,a,12x,a)') 'E (Ha)', 'Error (Ha)'

   ! +-----------------------------------------------------------------------+
   ! |                    READ INTEGRALS & DENSITIES                         |
   ! +-----------------------------------------------------------------------+
   call section(ouf, 'Reading integrals file & SIRIFC', ' ', termwidth)

   ! Read the integrals
   call qopen(intfn, 'formatted', 'sequential', intlu)
   call read_integrals(intlu, norb, repnuc, H1MO, H2MO)
   call qclose(intlu, intfn)

   ! Read the SIRIFC file
   call qopen(sirifn, 'unformatted', 'sequential', sirilu)
   call read_sirifc
   call qclose(sirilu, sirifn)

   if (tstdbg) then 
      call print_matrix('H1MO', H1MO(1:pmatsize,1:pmatsize), ouf, dfmt)
      call print_matrix('D1MO', D1MO(1:pmatsize,1:pmatsize), ouf, dfmt)
      call print_matrix('H2MO', H2MO(1,1,1:pmatsize,1:pmatsize), ouf, dfmt)
      call print_matrix('D2MO', D2MO(1,1,1:pmatsize,1:pmatsize), ouf, dfmt)
   end if

   ! Check that D1MO is read correctly
   call check_virtorbs_D1(D1MO, nocct, norbt, ouf)

   !
   ! Extract occupation numbers
   !
   call is_diag(D1MO, fthres, diag)
   if (diag) then
      write(ouf,'(a,/)') 'D1MO is already diagonal. Extracting ONs...'
      call extract_diag(D1MO, NONMO)
   else
      write(ouf,'(a,/)') 'D1MO is not diagonal. Diagonalizing...'
      call diag_2D_mat(D1MO, .true., NONMO)
      write(ouf,'(a,/)') 'Extracting ONs...'
      call extract_diag(D1MO, NONMO)
   end if
   if (tstdbg) then 
      write(ouf,'(a)') 'Occupied MO orbitals:'
      write(ouf,'(a)') '#  i  n_i'
      j = 1
      do i = 1, size(NONMO)
         if (NONMO(i).gt.fthres) then
            write(ouf,'(2(i0,2x),f4.2)') j, i, NONMO(i)
            j = j + 1
         end if
      end do
      write(ouf,*)
   end if

   ! +-----------------------------------------------------------------------+
   ! |                  COMPUTE EXACT ENERGY (NONSYMMETRIC)                  |
   ! +-----------------------------------------------------------------------+

   ! Calculate CASSCF energy as in Dalton to check integrals
   call section(ouf, 'Calculating CASSCF (nonsymmetric)', ' ', termwidth)
   call calc_casscf_dalton

   ! Calculate CASSCF energy from non-symmetric densities
   write(ouf,'(a)') 'Energy calculated with gamma, Gamma (nonsymmetric):'
   call calc_Eoe(D1MO, H1MO, nisht, nasht, Eoe_ina, Eoe_act, Eoe)
   call calc_Eee(D2MO, H2MO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   Eoe_ina_exact = Eoe_ina 
   Eoe_act_exact = Eoe_act 
   Eee_ina_exact = Eee_ina 
   Eee_act_exact = Eee_act 
   Eee_cross_exact = Eee_cross 
   Eoe_exact = Eoe
   Eee_exact = Eee
   call print_energy
   call output_energy(eolu, 'SIRI', siri_emcscf, 0.d0)
   call output_energy(eolu, 'D1D2', E_MCSCF, dabs(siri_emcscf-E_MCSCF))
   !
   if (tstdbg) then
      call calc_casscf_nonsymmetric
   end if

   ! +-----------------------------------------------------------------------+
   ! |                          SYMMETRIZE MATRICES                          |
   ! +-----------------------------------------------------------------------+
   call section(ouf, 'Symmetrizing 2-RDM', ' ', termwidth)
   D2MOsym = D2MO
   call symmetrize_D2(D2MOsym, nisht, nocct, tstdbg, ouf)
   call check_virtorbs_D2(D2MOsym, nocct, norbt, ouf)
   if (tstdbg) then 
      write(ouf,*)
      call print_matrix('D2MO', D2MO(1,1,1:pmatsize,1:pmatsize), ouf, dfmt)
      call print_matrix('D2MOsym',D2MOsym(1,1,1:pmatsize,1:pmatsize),ouf,dfmt)
   end if
   ! TODO
   ! D2MO = D2MOsym

   nelec = nactel + 2*nisht
   call checknorm_D1(D1MO, nelec, fthres, ouf)
   call checknorm_D2(D2MOsym, nelec, fthres, ouf)
   
   ! +-----------------------------------------------------------------------+
   ! |                  COMPUTE EXACT ENERGY (SYMMETRIC)                     |
   ! +-----------------------------------------------------------------------+
   call section(ouf, 'Calculating CASSCF (symmetric)', ' ', termwidth)

   write(ouf,'(a)') 'Energy calculated with gamma, Gamma (symmetrized):'
   call calc_Eee_sym(D2MOsym, H2MO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy

   ! +-----------------------------------------------------------------------+
   ! |               BBCn APPROXIMATION in SPACIAL ORBITALS                  |
   ! +-----------------------------------------------------------------------+
   call section(ouf, '(Müller) Buijse-Baerends corrected (BBCn)', ' ', termwidth)

   !
   ! Using spacial orbitals
   !
   write(ouf,'(a)') '--Spacial orbitals'

   write(ouf,'(a)') 'BBC1 :'
   call Eee_BBC_spacial('BBC1', NONMO/2.d0, H2MO, nisht, nasht, &
   & Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   call output_energy(eolu, 'BBC1spac', E_MCSCF, dabs(siri_emcscf-E_MCSCF))

   write(ouf,'(a)') 'BBC2 :'
   call Eee_BBC_spacial('BBC2', NONMO/2.d0, H2MO, nisht, nasht, &
   & Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   call output_energy(eolu, 'BBC2spac', E_MCSCF, dabs(siri_emcscf-E_MCSCF))

   write(ouf,'(a)') 'BBC3 :'
   call Eee_BBC_spacial('BBC3', NONMO/2.d0, H2MO, nisht, nasht, &
   & Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   call output_energy(eolu, 'BBC3spac', E_MCSCF, dabs(siri_emcscf-E_MCSCF))

   write(ouf,'(a)') 'BBC3M:'
   call Eee_BBC_spacial('BBC3M', NONMO/2.d0, H2MO, nisht, nasht, &
   & Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   call output_energy(eolu, 'BBC3Mspac', E_MCSCF, dabs(siri_emcscf-E_MCSCF))
   write(ouf,*)

   ! +-----------------------------------------------------------------------+
   ! |                        PRUEBA                                      |
   ! +-----------------------------------------------------------------------+
   call section(ouf, 'PRUEBA', ' ', termwidth)
   ! write(ouf,*) (2.d0 * NONMO(1)**2 - NONMO(1)) * H2MO(1,1,1,1)

   write(ouf,'(a)') 'LS'
   call calc_Eee_LS(NONMO, H2MO, nisht, nasht, &
   & Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   write(ouf,'(a)') 'BBC1'
   call calc_Eee_BBC('BBC1', NONMO, H2MO, nisht, nasht, &
   & Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   call output_energy(eolu, 'BBC1', E_MCSCF, dabs(siri_emcscf-E_MCSCF))
   write(ouf,'(a)') 'BBC2'
   call calc_Eee_BBC('BBC2', NONMO, H2MO, nisht, nasht, &
   & Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   call output_energy(eolu, 'BBC2', E_MCSCF, dabs(siri_emcscf-E_MCSCF))
   write(ouf,'(a)') 'BBC3'
   call calc_Eee_BBC('BBC3', NONMO, H2MO, nisht, nasht, &
   & Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   call output_energy(eolu, 'BBC3', E_MCSCF, dabs(siri_emcscf-E_MCSCF))
   write(ouf,'(a)') 'BBC3M'
   call calc_Eee_BBC('BBC3M', NONMO, H2MO, nisht, nasht, &
   & Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   call output_energy(eolu, 'BBC3M', E_MCSCF, dabs(siri_emcscf-E_MCSCF))

   ! +-----------------------------------------------------------------------+
   ! |                        APPROXIMATED 2RDMS                             |
   ! +-----------------------------------------------------------------------+
   call section(ouf, 'Approximated 2-RDms', ' ', termwidth)
   write(ouf,'(a)') '>> Approximated 2-RDms and normalization check for:'
   !
   ! LS
   !
   write(ouf,'(a)') '-- LS:'
   call D2_LS(NONMO, D2LSMO)
   call checknorm_D2(D2LSMO, nelec, fthres, ouf)
   !
   ! BBC1
   !
   write(ouf,'(a)') '-- BBC1:'
   call D2_BBC('BBC1', NONMO, D2BBC1MO)
   call checknorm_D2(D2BBC1MO, nelec, fthres, ouf)
   !
   ! BBC2
   !
   write(ouf,'(a)') '-- BBC2:'
   call D2_BBC('BBC2', NONMO, D2BBC2MO)
   call checknorm_D2(D2BBC2MO, nelec, fthres, ouf)
   !
   ! BBC3
   !
   write(ouf,'(a)') '-- BBC3:'
   call D2_BBC('BBC3', NONMO, D2BBC3MO)
   call checknorm_D2(D2BBC3MO, nelec, fthres, ouf)
   !
   ! BBC3
   !
   write(ouf,'(a)') '-- BBC3M:'
   call D2_BBC('BBC3M', NONMO, D2BBC3MMO)
   call checknorm_D2(D2BBC3MMO, nelec, fthres, ouf)

   if (tstdbg) then 
      call print_matrix('D2MO',      D2MO(1,1,1:pmatsize,1:pmatsize), ouf, dfmt)
      call print_matrix('D2BBC1MO',  D2BBC1MO(1,1,1:pmatsize,1:pmatsize), ouf, dfmt)
      call print_matrix('D2BBC2MO',  D2BBC2MO(1,1,1:pmatsize,1:pmatsize), ouf, dfmt)
      call print_matrix('D2BBC3MO',  D2BBC3MO(1,1,1:pmatsize,1:pmatsize), ouf, dfmt)
      call print_matrix('D2BBC3MMO', D2BBC3MMO(1,1,1:pmatsize,1:pmatsize), ouf, dfmt)
   end if

   ! +-----------------------------------------------------------------------+
   ! |                          MINKOWSKI DISTANCE                            |
   ! +-----------------------------------------------------------------------+
   call section(ouf, 'Minkowski distance for the approximated 2-RDMs', ' ', termwidth)

   write(ouf,'(a,i0,":")') '>> Minkowski metric (distance) of order ',minkord

   if (tstdbg) then 
      call Minkowski_distance_4D(minkord, D2MO, D2MO, mindis)
      write(ouf,'(a, f8.4)') 'Exact   :', mindis
   end if

   call Minkowski_distance_4D(minkord, D2MO, D2LSMO, mindis_LS)
   write(ouf,'(a, f8.4)') 'LS   :', mindis_LS

   call Minkowski_distance_4D(minkord, D2MO, D2BBC1MO, mindis_BBC1)
   write(ouf,'(a, f8.4)') 'BBC1 :', mindis_BBC1

   call Minkowski_distance_4D(minkord, D2MO, D2BBC2MO, mindis_BBC2)
   write(ouf,'(a, f8.4)') 'BBC2 :', mindis_BBC2

   call Minkowski_distance_4D(minkord, D2MO, D2BBC3MO, mindis_BBC3)
   write(ouf,'(a, f8.4)') 'BBC3 :', mindis_BBC3

   call Minkowski_distance_4D(minkord, D2MO, D2BBC3MMO, mindis_BBC3M)
   write(ouf,'(a, f8.4)') 'BBC3M:', mindis_BBC3M
   write(ouf,*)

   ! +-----------------------------------------------------------------------+
   ! |                      ENERGIES FROM APPROX 2RDMS                       |
   ! +-----------------------------------------------------------------------+
   call section(ouf, 'Energies with approximated 2-RDMs', ' ', termwidth)
   write(ouf,'(a)') '>> Energy from the approximated 2-RDM for:'
   !
   ! LS
   !
   write(ouf,'(a)') '-- LS:'
   ! D2LSMO(1,1,1,1) = 5.d0
   ! D2LSMO(3,3,3,3) = 5.d0
   ! D2LSMO(4,4,4,4) = 5.d0
   ! D2LSMO(5,5,5,5) = 5.d0
   ! D2LSMO(6,6,6,6) = 5.d0
   call calc_Eee(D2LSMO, H2MO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   E_LS = Eee + Eoe + repnuc
   call print_energy
   call output_energy(eolu, 'LS D2', E_MCSCF, dabs(siri_emcscf-E_MCSCF))
   !
   ! BBC1
   !
   write(ouf,'(a)') '-- BBC1:'
   call calc_Eee(D2BBC1MO, H2MO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   E_BBC1 = Eee + Eoe + repnuc
   call print_energy
   call output_energy(eolu, 'BBC1 D2', E_MCSCF, dabs(siri_emcscf-E_MCSCF))
   !
   ! BBC2
   !
   write(ouf,'(a)') '-- BBC2:'
   call calc_Eee(D2BBC2MO, H2MO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   E_BBC2 = Eee + Eoe + repnuc
   call print_energy
   call output_energy(eolu, 'BBC2 D2', E_MCSCF, dabs(siri_emcscf-E_MCSCF))
   !
   ! BBC3
   !
   write(ouf,'(a)') '-- BBC3:'
   call calc_Eee(D2BBC3MO, H2MO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   E_BBC3 = Eee + Eoe + repnuc
   call print_energy
   call output_energy(eolu, 'BBC3 D2', E_MCSCF, dabs(siri_emcscf-E_MCSCF))
   !
   ! BBC3
   !
   write(ouf,'(a)') '-- BBC3M:'
   call calc_Eee(D2BBC3MMO, H2MO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   E_BBC3M = Eee + Eoe + repnuc
   call print_energy
   call output_energy(eolu, 'BBC3M D2', E_MCSCF, dabs(siri_emcscf-E_MCSCF))

   ! +-----------------------------------------------------------------------+
   ! |                              DUMP RESULTS                             |
   ! +-----------------------------------------------------------------------+
   call qopen('plot_data', 'formatted', 'sequential', epolu)

   write(epolu,'(a)')&
   &'            LS            BBC1          BBC2          BBC3          BBC3M'
   write(epolu,'(a,*(f14.6))')&
   &'E (Ha) :', E_LS, E_BBC1, E_BBC2, E_BBC3, E_BBC3M
   write(epolu,'(a,*(d14.2))')&
   &'ΔE (Ha):', &
   &dabs(E_LS-siri_emcscf),&
   &dabs(E_BBC1-siri_emcscf),&
   &dabs(E_BBC2-siri_emcscf),&
   &dabs(E_BBC3-siri_emcscf),&
   &dabs(E_BBC3M-siri_emcscf)
   write(epolu,'(a,*(f14.6))')&
   &'Mink. d:', mindis_LS, mindis_BBC1, mindis_BBC2, mindis_BBC3, mindis_BBC3M
   
   call qclose(epolu, 'plot_data')

   ! +-----------------------------------------------------------------------+
   ! |                     TRANSFORM TO SPIN-ORBITALS                        |
   ! +-----------------------------------------------------------------------+
   ! call section(ouf, 'Change to spin-orbital basis', ' ', termwidth)
   !
   ! ! Transform D^(1), D^(2), H1, H2 from MO to spin-orbital basis
   ! call MO_to_SO_D1(D1MO, D1SO)
   ! call MO_to_SO_D2(D2MO, D2SO)  ! TODO: usar D2 simetrizada
   ! call MO_to_SO_H1(H1MO, H1SO)
   ! call MO_to_SO_H2(H2MO, H2SO)
   !
   ! if (tstdbg) then 
   !    call print_matrix('H1SO', H1SO(1:2*pmatsize,1:2*pmatsize), ouf, dfmt)
   !    call print_matrix('D1SO', D1SO(1:2*pmatsize,1:2*pmatsize), ouf, dfmt)
   !    call print_matrix('H2SO', H2SO(1,1,1:2*pmatsize,1:2*pmatsize), ouf, dfmt)
   !    call print_matrix('D2SO', D2SO(1,1,1:2*pmatsize,1:2*pmatsize), ouf, dfmt)
   ! end if
   !
   ! ! Check normalization
   ! call checknorm_D1(D1SO, nelec, fthres, ouf)
   ! call checknorm_D2(D2SO, nelec, fthres, ouf)
   !
   ! ! Check SO occupations
   ! call is_diag(D1SO, fthres, diag)
   ! if (diag) then
   !    write(ouf,'(a,/)') 'D1SO is already diagonal. Extracting ONs...'
   !    call extract_diag(D1SO, NONSO)
   ! else
   !    write(ouf,'(a,/)') 'D1SO is not diagonal. Diagonalizing...'
   !    call diag_2D_mat(D1SO, .true., NONSO)
   !    write(ouf,'(a,/)') 'Extracting ONs...'
   !    call extract_diag(D1SO, NONSO)
   ! end if
   ! if (tstdbg) then 
   !    write(ouf,'(a)') 'Occupied SO orbitals:'
   !    write(ouf,'(a)') '#  i  n_i'
   !    j = 1
   !    do i = 1, size(NONSO)
   !       if (NONSO(i).gt.fthres) then
   !          write(ouf,'(2(i0,2x),f4.2)') j, i, NONSO(i)
   !          j = j + 1
   !       end if
   !    end do
   !    write(ouf,*)
   ! end if
   !
   ! !
   ! ! Recalculate the energy to check that the conversion is right
   ! !
   ! ! TODO
   ! if (tstdbg) then
   !    call calc_Eee_SO(D2SO, H2SO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   !    call print_energy
   ! end if

   ! +-----------------------------------------------------------------------+
   ! |                  BBCn APPROXIMATION in SPIN-ORBITALS                  |
   ! +-----------------------------------------------------------------------+
   !
   ! Using spin orbitals
   !
   ! write(ouf,'(a)') 'Spin orbitals'
   !
   ! call calc_Eee_BBC('BBC1', NONSO, H2SO, nisht, nasht, &
   ! & Eee_ina, Eee_act, Eee_cross, Eee)
   ! write(ouf,'(a)') 'BBC1:' ; call print_energy
   ! !
   ! call calc_Eee_BBC('BBC2', NONSO, H2SO, nisht, nasht, &
   ! & Eee_ina, Eee_act, Eee_cross, Eee)
   ! write(ouf,'(a)') 'BBC2:' ; call print_energy
   ! !
   ! call calc_Eee_BBC('BBC3', NONSO, H2SO, nisht, nasht, &
   ! & Eee_ina, Eee_act, Eee_cross, Eee)
   ! write(ouf,'(a)') 'BBC3:' ; call print_energy
   ! !
   ! call calc_Eee_BBC('BBC3M', NONSO, H2SO, nisht, nasht, &
   ! & Eee_ina, Eee_act, Eee_cross, Eee)
   ! write(ouf,'(a)') 'BBC3M:' ; call print_energy
   ! !
   ! !
   ! ! OLD SUBROUTINES
   ! !
   ! call Eee_BBC('BBC1', NONSO, H2SO, Eee)
   ! !
   ! E_elec = Eoe + Eee
   ! E_MCSCF = E_elec + repnuc
   ! write(ouf,'(a,2f16.10)') 'BBC1 MCSCF energy, error    :', E_MCSCF, &
   ! & dabs(siri_emcscf-E_MCSCF)
   ! call flush(ouf)
   !
   ! call Eee_BBC('BBC2', NONSO, H2SO, Eee)
   ! !
   ! E_elec = Eoe + Eee
   ! E_MCSCF = E_elec + repnuc
   ! write(ouf,'(a,2f16.10)') 'BBC2 MCSCF energy, error    :', E_MCSCF, &
   ! & dabs(siri_emcscf-E_MCSCF)
   !
   ! call Eee_BBC('BBC3', NONSO, H2SO, Eee)
   ! !
   ! E_elec = Eoe + Eee
   ! E_MCSCF = E_elec + repnuc
   ! write(ouf,'(a,2f16.10)') 'BBC3 MCSCF energy, error    :', E_MCSCF, &
   ! & dabs(siri_emcscf-E_MCSCF)
   !
   ! call Eee_BBC('BBC3M', NONSO, H2SO, Eee)
   ! !
   ! E_elec = Eoe + Eee
   ! E_MCSCF = E_elec + repnuc
   ! write(ouf,'(a,2f16.10)') 'BBC3M MCSCF energy, error    :', E_MCSCF, &
   ! & dabs(siri_emcscf-E_MCSCF)
   !
!
!
!    ! +-----------------------------------------------------------------------+
!    ! |                            LS 2-RDM MATRIX                            |
!    ! +-----------------------------------------------------------------------+
!    call section(ouf, 'LS approximated 2-RDMs', ' ', termwidth)
!
!    call D2_LS(NONMO, D2LSMO)
!    call D2_LS(NONSO, D2LSSO)
!
!    write(ouf,'(a)') 'Energy calculated LS D^(2):'
!    write(ouf,'(a)') 'D2LSMO:'
!    call calc_Eee(D2LSMO, H2MO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
!    call print_energy
!    write(ouf,'(a)') 'D2LSSO:'
!    call calc_Eee_SO(D2LSSO, H2SO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
!    call print_energy
!
!    ! +-----------------------------------------------------------------------+
!    ! |                            BBCn 2-RDM MATRIX                          |
!    ! +-----------------------------------------------------------------------+
!    call section(ouf, 'BBCn approximated 2-RDMs', ' ', termwidth)
!
!    call D2_BBC('BBC1', NONMO, D2BBC1MO)
!    call D2_BBC('BBC2', NONMO, D2BBC2MO)
!    ! call D2_BBC2(NONMO, D2BBC2MO)  ! TODO
!    call D2_BBC('BBC3', NONMO, D2BBC3MO)
!    call D2_BBC('BBC3M', NONMO, D2BBC3MMO)
!    !
!    call D2_BBC('BBC1', NONSO, D2BBC1SO)
!    call D2_BBC('BBC2', NONSO, D2BBC2SO)
!    call D2_BBC('BBC3', NONSO, D2BBC3SO)
!    call D2_BBC('BBC3M', NONSO, D2BBC3MSO)
!    !
!    ! Energies
!    !
!    write(ouf,'(a)') '--Energy calculated with BBCn D^(2):'
!    write(ouf,'(a)') 'BBC1SO:'
!    call calc_Eee(D2BBC1SO, H2SO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
!    call print_energy
!
!    write(ouf,'(a)') 'BBC1MO:'
!    call calc_Eee(D2BBC1MO, H2MO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
!    call print_energy
!    !
!    write(ouf,'(a)') 'BBC2:'
!    call calc_Eee(D2BBC2MO, H2MO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
!    call print_energy
!    !
!    write(ouf,'(a)') 'BBC3:'
!    call calc_Eee(D2BBC3MO, H2MO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
!    call print_energy
!    !
!    write(ouf,'(a)') 'BBC3M:'
!    call calc_Eee(D2BBC3MMO, H2MO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
!    call print_energy
!
!
!    ! +-----------------------------------------------------------------------+
!    ! |                          LS APPROXIMATION                             |
!    ! +-----------------------------------------------------------------------+
!    ! call section(ouf, 'Löwdin-Shull (LS)', ' ', termwidth)
!
!    ! Calculate energy
!    ! call  calc_Eee_LS(NONSO, H2SO, nisht, nasht, &
!    ! &                 Eee_ina, Eee_act, Eee_cross, Eee)
!    ! call print_energy
!    !
!    ! ! Calculate energy
!    ! call Eee_LS(NONSO, H2SO, Eee)
!    ! E_elec = Eoe + Eee
!    ! E_MCSCF = E_elec + repnuc
!    ! !
!    ! write(ouf,'(a)') 'LS Energy calculated with gamma:'
!    ! write(ouf,'(a,2f16.10,/)') 'MCSCF energy, error    :', E_MCSCF, &
!    ! & dabs(siri_emcscf-E_MCSCF)
!    ! call flush(ouf)
!    
   call qclose(eolu, 'ENERGIES')
   !
   stop
end program RDMFT
