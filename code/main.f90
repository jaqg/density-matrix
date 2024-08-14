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
      ! TODO: transform integrals to the same basis
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

   ! Store exact results
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

   ! Symmetrize & check the 2-RDM
   D2MOsym = D2MO
   call symmetrize_D2(D2MOsym, nisht, nocct, tstdbg, ouf)
   call check_virtorbs_D2(D2MOsym, nocct, norbt, ouf)

   if (tstdbg) then 
      write(ouf,*)
      call print_matrix('D2MO', D2MO(1,1,1:pmatsize,1:pmatsize), ouf, dfmt)
      call print_matrix('D2MOsym',D2MOsym(1,1,1:pmatsize,1:pmatsize),ouf,dfmt)
   end if

   ! Calculate the number of electrons
   nelec = nactel + 2*nisht

   ! Check the normalization conditions
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
   ! |                         MINKOWSKI DISTANCE                            |
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
   call section(ouf, 'Change to spin-orbital basis', ' ', termwidth)
   !
   ! Transform D^(1), D^(2), H1, H2 from MO to spin-orbital basis
   !
   !
   call MO_to_SO_D1(D1MO, D1SO)
   call MO_to_SO_D2(D2MO, D2SO)  ! TODO: usar D2 simetrizada
   call MO_to_SO_H1(H1MO, H1SO)
   call MO_to_SO_H2(H2MO, H2SO)
   ! test = 0.d0
   ! do i=1, norb
   !    test = test + D1MO(i,i) * H1MO(i,i)
   ! end do
   ! write(*,*) 'test MO: ', test
   ! test = 0.d0
   ! do i=1, 2*norb
   !    test = test + D1SO(i,i) * H1SO(i,i)
   ! end do
   ! write(*,*) 'test SO: ', test
   ! test = 0.d0
   ! do i=1, norb
   !    do j=1, norb
   !       test = test + D2MO(i,j,i,j) * H2MO(i,j,i,j)
   !    end do
   ! end do
   ! write(*,*) 'test MO: ', test
   ! test = 0.d0
   ! do i=1, norb
   !    do j=1, norb
   !       test = test + D2SO(i,j,i,j) * H2SO(i,j,i,j)
   !    end do
   ! end do
   ! write(*,*) 'test SO: ', test
   !
   if (tstdbg) then 
      call print_matrix('CMO',  CMO(1:pmatsize,1:pmatsize), ouf, dfmt)
      call print_matrix('H1SO', H1SO(1:2*pmatsize,1:2*pmatsize), ouf, dfmt)
      call print_matrix('D1SO', D1SO(1:2*pmatsize,1:2*pmatsize), ouf, dfmt)
      call print_matrix('H2SO', H2SO(1,1,1:2*pmatsize,1:2*pmatsize), ouf, dfmt)
      call print_matrix('D2SO', D2SO(1,1,1:2*pmatsize,1:2*pmatsize), ouf, dfmt)
   end if
   !
   ! Check normalization
   call checknorm_D1(D1SO, nelec, fthres, ouf)
   call checknorm_D2(D2SO, nelec, fthres, ouf)
   !
   ! Check SO occupations
   call is_diag(D1SO, fthres, diag)
   if (diag) then
      write(ouf,'(a,/)') 'D1SO is already diagonal. Extracting ONs...'
      call extract_diag(D1SO, NONSO)
   else
      write(ouf,'(a,/)') 'D1SO is not diagonal. Diagonalizing...'
      call diag_2D_mat(D1SO, .true., NONSO)
      write(ouf,'(a,/)') 'Extracting ONs...'
      call extract_diag(D1SO, NONSO)
   end if
   if (tstdbg) then 
      write(ouf,'(a)') 'Occupied SO orbitals:'
      write(ouf,'(a)') '#  i  n_i'
      j = 1
      do i = 1, size(NONSO)
         if (NONSO(i).gt.fthres) then
            write(ouf,'(2(i0,2x),f4.2)') j, i, NONSO(i)
            j = j + 1
         end if
      end do
      write(ouf,*)
   end if

   !
   ! Recalculate the energy to check that the conversion is right
   !
   ! call calc_Eee_SO(D2SO, H2SO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   ! call print_energy

   ! +-----------------------------------------------------------------------+
   ! |                        ENERGY IN SO BASIS                             |
   ! +-----------------------------------------------------------------------+
   call section(ouf, 'Energy in SO basis', ' ', termwidth)
   write(ouf,'(a)') '>> Energy functionals in SO basis for:'
   !
   ! LS
   !
   write(ouf,'(a)') '-- LS:'
   call calc_Eee_LS(NONSO, H2SO, 2*nisht, 2*nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   !
   ! BBC1
   !
   write(ouf,'(a)') '-- BBC1:'
   call calc_Eee_BBC('BBC1', NONSO, H2SO, 2*nisht, 2*nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   !
   ! BBC2
   !
   write(ouf,'(a)') '-- BBC2:'
   call calc_Eee_BBC('BBC2', NONSO, H2SO, 2*nisht, 2*nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   !
   ! BBC3
   !
   write(ouf,'(a)') '-- BBC3:'
   call calc_Eee_BBC('BBC3', NONSO, H2SO, 2*nisht, 2*nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   !
   ! BBC3M
   !
   write(ouf,'(a)') '-- BBC3M:'
   call calc_Eee_BBC('BBC3M', NONSO, H2SO, 2*nisht, 2*nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy

   ! +-----------------------------------------------------------------------+
   ! |                        APPROXIMATED 2RDMS                             |
   ! +-----------------------------------------------------------------------+
   call section(ouf, 'SO Approximated 2-RDms', ' ', termwidth)
   write(ouf,'(a)') '>> Approximated 2-RDms and normalization check for:'
   !
   ! LS
   !
   write(ouf,'(a)') '-- LS:'
   call D2_LS(NONSO, D2LSSO)
   call checknorm_D2(D2LSSO, nelec, fthres, ouf)
   !
   ! BBC1
   !
   write(ouf,'(a)') '-- BBC1:'
   call D2_BBC('BBC1', NONSO, D2BBC1SO)
   call checknorm_D2(D2BBC1SO, nelec, fthres, ouf)
   !
   ! BBC2
   !
   write(ouf,'(a)') '-- BBC2:'
   call D2_BBC('BBC2', NONSO, D2BBC2SO)
   call checknorm_D2(D2BBC2SO, nelec, fthres, ouf)
   !
   ! BBC3
   !
   write(ouf,'(a)') '-- BBC3:'
   call D2_BBC('BBC3', NONSO, D2BBC3SO)
   call checknorm_D2(D2BBC3SO, nelec, fthres, ouf)
   !
   ! BBC3
   !
   write(ouf,'(a)') '-- BBC3M:'
   call D2_BBC('BBC3M', NONSO, D2BBC3MSO)
   call checknorm_D2(D2BBC3MSO, nelec, fthres, ouf)

   if (tstdbg) then 
      call print_matrix('D2SO',      D2SO(1,1,1:pmatsize,1:pmatsize), ouf, dfmt)
      call print_matrix('D2BBC1SO',  D2BBC1SO(1,1,1:pmatsize,1:pmatsize), ouf, dfmt)
      call print_matrix('D2BBC2SO',  D2BBC2SO(1,1,1:pmatsize,1:pmatsize), ouf, dfmt)
      call print_matrix('D2BBC3SO',  D2BBC3SO(1,1,1:pmatsize,1:pmatsize), ouf, dfmt)
      call print_matrix('D2BBC3MSO', D2BBC3MSO(1,1,1:pmatsize,1:pmatsize), ouf, dfmt)
   end if

   ! +-----------------------------------------------------------------------+
   ! |                         MINKOWSKI DISTANCE                            |
   ! +-----------------------------------------------------------------------+
   call section(ouf, 'Minkowski distance for the approximated SO 2-RDMs', ' ', termwidth)

   write(ouf,'(a,i0,":")') '>> Minkowski metric (distance) of order ',minkord

   if (tstdbg) then 
      call Minkowski_distance_4D(minkord, D2SO, D2SO, mindis)
      write(ouf,'(a, f8.4)') 'Exact   :', mindis
   end if

   call Minkowski_distance_4D(minkord, D2SO, D2LSSO, mindis_LS)
   write(ouf,'(a, f8.4)') 'LS   :', mindis_LS

   call Minkowski_distance_4D(minkord, D2SO, D2BBC1SO, mindis_BBC1)
   write(ouf,'(a, f8.4)') 'BBC1 :', mindis_BBC1

   call Minkowski_distance_4D(minkord, D2SO, D2BBC2SO, mindis_BBC2)
   write(ouf,'(a, f8.4)') 'BBC2 :', mindis_BBC2

   call Minkowski_distance_4D(minkord, D2SO, D2BBC3SO, mindis_BBC3)
   write(ouf,'(a, f8.4)') 'BBC3 :', mindis_BBC3

   call Minkowski_distance_4D(minkord, D2SO, D2BBC3MSO, mindis_BBC3M)
   write(ouf,'(a, f8.4)') 'BBC3M:', mindis_BBC3M
   write(ouf,*)

   ! +-----------------------------------------------------------------------+
   ! |                      ENERGIES FROM APPROX 2RDMS                       |
   ! +-----------------------------------------------------------------------+
   call section(ouf, 'Energies with approximated SO 2-RDMs', ' ', termwidth)
   write(ouf,'(a)') '>> Energy from the approximated 2-RDM for:'
   !
   ! LS
   !
   write(ouf,'(a)') '-- LS:'
   call calc_Eee(D2LSSO, H2SO, 2*nisht, 2*nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   !
   ! BBC1
   !
   write(ouf,'(a)') '-- BBC1:'
   call calc_Eee(D2BBC1SO, H2SO, 2*nisht, 2*nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   !
   ! BBC2
   !
   write(ouf,'(a)') '-- BBC2:'
   call calc_Eee(D2BBC2SO, H2SO, 2*nisht, 2*nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   !
   ! BBC3
   !
   write(ouf,'(a)') '-- BBC3:'
   call calc_Eee(D2BBC3SO, H2SO, 2*nisht, 2*nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy
   !
   ! BBC3
   !
   write(ouf,'(a)') '-- BBC3M:'
   call calc_Eee(D2BBC3MSO, H2SO, 2*nisht, 2*nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   call print_energy

   ! +-----------------------------------------------------------------------+
   ! |                        OPTIMIZATION ALGORITHM                         |
   ! +-----------------------------------------------------------------------+

   ! +-----------------------------------------------------------------------+
   ! |                                   END                                 |
   ! +-----------------------------------------------------------------------+
   call qclose(eolu, 'ENERGIES')
   !
   stop
end program RDMFT
