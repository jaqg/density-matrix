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

   ! +-----------------------------------------------------------------------+
   ! |                    READ INTEGRALS & DENSITIES                         |
   ! +-----------------------------------------------------------------------+
   call section(ouf, 'Reading integrals file & SIRIFC', ' ', termwidth)

   ! Read the integrals
   call qopen(intfn, 'formatted', 'sequential', intlu)
   call read_integrals(intlu, norb, repnuc, H1MO, H2MO)
   call qclose(intlu, intfn)


   ! call read_integrals_alfredo

   ! Read the SIRIFC file
   call qopen(sirifn, 'unformatted', 'sequential', sirilu)
   call read_sirifc
   call qclose(sirilu, sirifn)

   ! TODO: check normalizacion de las RDMS 
   nelec = nactel + 2*nisht
   write(*,*) 'nactel =', nactel
   write(*,*) 'nelec =', nelec
   ! call checknorm_D1(D1MO, nelec, fthres, ouf)

   ! TODO: simetrizar 2RDM (?)
   
   if (tstdbg) then 
      call print_matrix('H1MO', H1MO(1:pmatsize,1:pmatsize), ouf, dfmt)
      call print_matrix('D1MO', D1MO(1:pmatsize,1:pmatsize), ouf, dfmt)
   end if

   ! +-----------------------------------------------------------------------+
   ! |                       COMPUTE EXACT ENERGY                            |
   ! +-----------------------------------------------------------------------+

   ! Calculate CASSCF energy as in Dalton to check integrals
   call section(ouf, 'Calculating CASSCF', ' ', termwidth)
   call calc_casscf_dalton

   ! Calculate CASSCF energy from non-symmetric densities
   call calc_Eoe(D1MO, H1MO, nisht, nasht, Eoe_ina, Eoe_act, Eoe)
   call calc_Eee(D2MO, H2MO, nisht, nasht, Eee_ina, Eee_act, Eee_cross, Eee)
   E_ina = Eoe_ina + Eee_ina
   E_act = Eoe_act + Eee_cross + Eee_act
   E_elec = Eoe + Eee
   E_MCSCF = E_elec + repnuc

   write(ouf,'(a)') 'Energy calculated with gamma, Gamma'
   write(ouf,'(a,2f16.10)') 'Inactive energy, error :', E_ina, &
   & dabs(E_ina-siri_einactiv)
   write(ouf,'(a,2f16.10)') 'act-ina contribution   :', Eoe_act + Eee_cross
   write(ouf,'(a,2f16.10)') 'act-act contribution   :', Eee_act 
   write(ouf,'(a,2f16.10)') 'Active energy, error   :', E_act, &
   & dabs(siri_eactiv-E_act)
   write(ouf,'(a,2f16.10,/)') 'MCSCF energy, error   :', E_MCSCF, &
   & dabs(siri_emcscf-E_MCSCF)
   call flush(ouf)
   !
   call calc_casscf_nonsymmetric

   !
   stop
end program RDMFT
