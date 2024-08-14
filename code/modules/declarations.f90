! +--------------------------------------------+
! | Author: Jose Antonio Quinonero Gris        |
! | Creation date: Tuesday 09:11:31 16/01/2024 |
! +--------------------------------------------+
module declarations
   ! 
   ! Module to declare variables
   !
   implicit none
   !
   logical :: tstdbg, diag
   !
   character(len=8) :: stars(3), label
   character(len=80) :: rfmt, dfmt, intfn, sirifn
   !
   integer(kind=8) :: i, j, k, l, n, m
   integer :: ierr, ouf, eolu, epolu, intlu, sirilu, termwidth, pmatsize, minkord 
   integer(kind=8) :: ij, ispin, istate, koff, lsym, ms2, n2orbt, &
   & nacorb, nactel, nasht, nbast, ncdets, ncmot, nconf, ni, ninorb, nisht, &
   & nj, nnashx, nnashy, nnorbt, nocct, norbt, nseorb, nsym, nu, numorb, &
   & nv, nwoph, nwopt, nx, ny, nelec
   integer(kind=8) :: p,q,r,s,u,v,x,y,uv,xy
   integer(kind=4) :: norb 
   !
   real(kind=8) :: fthres
   real(kind=8) :: eact, eact1, eact2, siri_eactiv, eina, siri_emcscf, siri_einactiv, erract, &
   & errina, repnuc, siri_repnuc, xint
   real(kind=8) :: E_act, E_ina, E_cross, E_elec, E_MCSCF 
   real(kind=8) :: Eoe_ina, Eee_ina, Eoe_act, Eee_act, Eee_cross, Eoe, Eee
   real(kind=8) :: Eoe_ina_exact, Eee_ina_exact, Eoe_act_exact, Eee_act_exact,&
   & Eee_cross_exact, E_cross_exact, Eoe_exact, Eee_exact
   real(kind=8) :: E_LS, E_BBC1, E_BBC2, E_BBC3, E_BBC3M
   real(kind=8) :: mindis, mindis_LS, mindis_BBC1, mindis_BBC2, mindis_BBC3, &
   & mindis_BBC3M
   !
   real(kind=8), allocatable, dimension (:,:) :: oneint, den1
   real(kind=8), allocatable, dimension (:,:,:,:) :: twoint, rdm2
   real(kind=8), allocatable, dimension (:,:) :: fockin
   real(kind=8), allocatable, dimension (:) :: qmat, scr_cmo, scr
   real(kind=8), allocatable, dimension (:,:) :: CMO, H1MO, D1MO, H1SO, D1SO
   real(kind=8), allocatable, dimension (:,:,:,:) :: H2MO, D2MO, D2MOsym, &
   & H2SO, D2SO
   real(kind=8), allocatable, dimension (:,:,:,:) :: D2LSMO, D2LSSO
   real(kind=8), allocatable, dimension (:,:,:,:) :: D2BBC1MO, D2BBC1SO, &
   & D2BBC2MO, D2BBC2SO, D2BBC3MO, D2BBC3SO, D2BBC3MMO, D2BBC3MSO
   real(kind=8), allocatable, dimension (:) :: NONMO, NONSO
   real(kind=8) :: test 
   !
   contains
       subroutine parameters
        implicit none

        ! Print debug information
        tstdbg = .false.

        ! floating-point numbers threshold to consider it 0, i.e.
        ! a -> 0 if a < fthres
        fthres = 1.d-10

        ! Width of the terminal/file to be written on
        termwidth = 60

        ! Logic unit of the output file
        ouf = 6

        ! D format
        dfmt = '(*(d10.2))'
        rfmt = '(*(f10.4))'

        ! Matrix size to print: print A(1:pmatsize,1:pmatsize)
        pmatsize = 4

        ! Minkowski metric order, p
        minkord = 4

        return
    end subroutine parameters 
 
end module declarations
