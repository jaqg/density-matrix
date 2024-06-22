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
   logical :: tstdbg
   !
   character(len=8) :: stars(3), label
   character(len=80) :: rfmt, dfmt, intfn, sirifn
   !
   integer(kind=8) :: i, j, k, l, n, m
   integer :: ierr, ouf, intlu, sirilu, termwidth, pmatsize, minkord 
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
   !
   real(kind=8), allocatable, dimension (:,:) :: oneint, den1, H1MO, D1MO
   real(kind=8), allocatable, dimension (:,:,:,:) :: twoint, rdm2, H2MO, D2MO
   real(kind=8), allocatable, dimension (:,:) :: fockin
   real(kind=8), allocatable, dimension (:) :: qmat, scr
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
