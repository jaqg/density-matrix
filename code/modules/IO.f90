! +--------------------------------------------+
! | Author: Jose Antonio Quinonero Gris        |
! | Creation date: Tuesday 17:27:38 30/04/2024 |
! +--------------------------------------------+
module IO
   ! 
   ! Input/Output module
   !
   use declarations
   use strings
   use matmod
   ! use density_matrices

   implicit none
   contains

   integer(kind=8) function indij(i,j)
      implicit none
      integer(kind=8) :: i, j

      indij = max(i,j)*(max(i,j)-1)/2 + min(i,j)

   end function indij


   logical function file_exist(filename)
      !
      ! Function to check if file 'filename' exists
      !
      implicit none
      character(len=*), intent(in) :: filename
      logical :: exists

      inquire(file=filename, exist=exists)
      file_exist = exists

   end function file_exist

   subroutine qopen(fn, rformat, raccess, lu)
      !
      ! Subroutine to open file with filename 'fn' in a new unit 'lu'
      !
      implicit none
      character(len=*), intent(in) :: fn, rformat, raccess
      integer, intent(out) :: lu
      !
      integer :: ios 

      open(newunit=lu, file=trim(fn), form=trim(rformat), access=trim(raccess),&
      & iostat=ios)
      if (ios /= 0) stop 'ERROR qopen: Error opening file ' // trim(fn)

      !
      return
   end subroutine qopen 

   subroutine qclose(lu, fn)
      !
      ! Subroutine to close the file 'fn' opened on the unit 'lu'
      !
      implicit none
      integer, intent(in) :: lu
      character(len=*), intent(in) :: fn 
      !
      integer :: ios 

      close(unit=lu, iostat=ios, status="keep")
      if (ios /= 0) stop "ERROR qclose: Error closing file " // trim(fn)

      !
   end subroutine qclose

   subroutine read_input(ifn, sirifn, intfn)
      !
      ! ifn: input file name
      ! sirifn: SIRIFC file name to read the density matrices from
      ! intfn: integrals file name
      !
      implicit none
      character(len=*), intent(in) :: ifn
      character(len=80), intent(out) :: sirifn, intfn
      !
      integer :: ios, nlu 
      logical :: file_exists
      character(len=80) :: datafolder 

      ! Open input file

      open(newunit=nlu, file=trim(ifn), iostat=ios)
      if (ios /= 0) stop "read_input: Error opening file trim(ifn)"

      ! Read input file
      read(nlu,*)       ! skip line
      read(nlu,*) datafolder
      read(nlu,*)
      read(nlu,*) sirifn
      read(nlu,*)
      read(nlu,*) intfn

      ! Add full path to the filenames
      sirifn = trim(datafolder) // '/' // trim(sirifn)
      intfn = trim(datafolder) // '/' // trim(intfn)

      ! Close file
      close(unit=nlu, iostat=ios, status="keep")
      if (ios /= 0) stop "read_input: Error closing file unit nlu"

      ! Check if files exist
      file_exists = file_exist(sirifn)
      if (.not.file_exists) then
         stop 'ERROR read_input: File ' // trim(sirifn) // ' does not exist'
      endif
      !
      file_exists = file_exist(intfn)
      if (.not.file_exists) then
         stop 'ERROR read_input: File ' // trim(intfn) // ' does not exist'
      endif

      return
   end subroutine read_input 

   subroutine read_integrals(intlu, norb, repnuc, H1, H2)
      !
      ! Subroutine to read the integrals
      !
      implicit none
      integer, intent(in) :: intlu
      integer, intent(out) :: norb
      real(kind=8), intent(out) :: repnuc 
      real(kind=8), dimension(:,:), allocatable, intent(out) :: H1
      real(kind=8), dimension(:,:,:,:), allocatable, intent(out) :: H2
      !
      logical :: found_norb
      integer :: i, j, k, l, pos, ios
      real(kind=8) :: xint 
      character(len=100) :: line 

      ! Rewind integrals file
      rewind(intlu)

      ! Read dimensions of H2
      found_norb = .false.

      ! Read the file
      read(intlu, *)  ! skip first line
      read(intlu, '(A)', iostat=ios) line
      if (ios /= 0) stop 'ERROR read_integrals: Error reading file.'

      ! Search for "NORB="
      pos = index(line, 'NORB=')
      if (pos /= 0) then
         read(line(pos+5:), *) norb  ! Read NORB
      end if

      ! Skip following 3 lines
      do i=1,3; read(intlu,*) ; end do

      ! Allocate H1
      allocate(H1(norb,norb), stat=ios)
      if (ios /= 0) stop "ERROR read_integrals: Error in allocation of H1"
      H1 = 0.d0

      ! Allocate H2
      allocate(H2(norb,norb,norb,norb), stat=ios)
      if (ios /= 0) stop "ERROR read_integrals: Error in allocation of H2"
      H2 = 0.d0

      ! Read the triangular part of H2mat
      i = -1
      j = -1
      k = -1
      l = -1
      dl1:do
         ! read(intlu,*,iostat=ios) xint, i, j, k, l
         ! if (ios /= 0) stop 'ERROR read_integrals: Error reading file.'
         read(intlu,*) xint, i, j, k, l

         if (l .ne. 0) then
            H2(i,j,k,l) = xint
            H2(i,j,l,k) = xint
            H2(j,i,k,l) = xint
            H2(j,i,l,k) = xint
            H2(k,l,i,j) = xint
            H2(k,l,j,i) = xint
            H2(l,k,i,j) = xint
            H2(l,k,j,i) = xint
         elseif (j .ne. 0) then
            H1(i,j) = xint
            H1(j,i) = xint
         else
            repnuc = xint
            exit dl1
         end if
      end do dl1

      !
      return
   end subroutine read_integrals 

   subroutine read_sirifc
      !
      ! Read SIRIFC file
      !
      rewind(sirilu)
      !
      read(sirilu) stars, label
      if ((stars(1) .ne. '********')&
         &    .or. (label(1:7) .ne. 'SIR IPH')) stop 'Error in SIRIFC'
      !
      read(sirilu) siri_repnuc, siri_einactiv, siri_eactiv, siri_emcscf, istate, ispin,&
      &         nactel, lsym, ms2
      !
      if (dabs(siri_repnuc-repnuc) .gt. 1.0d-8) then
         write(ouf,*) 'ERROR: Inconsistency between input and SIRIFC'
         write(ouf,*) 'Nuclear repulsion in input  :',repnuc
         write(ouf,*) 'Nuclear repulsion in SIRIFC :',siri_repnuc
         stop
      end if
      !
      write(ouf,'(a)') 'Data from SIRIFC:'
      write(ouf,'(a,f16.8)') 'Nuclear repulsion :', siri_repnuc
      write(ouf,'(a,f16.8)') 'Inactive energy   :', siri_einactiv
      write(ouf,'(a,f16.8)') 'Active energy     :', siri_eactiv
      write(ouf,'(a,f16.8)') 'MCSCF energy      :', siri_emcscf
      write(ouf,*)
      call flush(ouf)
      !
      read(sirilu) nisht, nasht, nocct, norbt, nbast, nconf,&
      &         nwopt, nwoph, ncdets, ncmot, nnashx, nnashy,&
      &         nnorbt, n2orbt, nsym
      !
      write(*,*) 'NSYM SIRIFC: ', nsym
      ! if (nsym .ne. 1) stop 'Symmetry not implemented yet'
      !
      if (norbt .ne. norb) then
         write(ouf,*) 'ERROR: Inconsistency between input and SIRIFC'
         write(ouf,*) 'Number of molecular orbitals in SIRIFC :', norbt
         write(ouf,*) 'Number of molecular orbitals in input  :', norb
         stop
      end if
      !
      ! if (nisht .ne. ninorb) then
      !    write(ouf,*) 'ERROR: Inconsistency between input and SIRIFC'
      !    write(ouf,*) 'Number of inactive orbitals in SIRIFC :', nisht
      !    write(ouf,*) 'Number of inactive orbitals in input  :', ninorb
      !    stop
      ! end if
      ! !
      ! if (nasht .ne. nacorb) then
      !    write(ouf,*) 'ERROR: Inconsistency between input and SIRIFC'
      !    write(ouf,*) 'Number of active orbitals in SIRIFC :', nasht
      !    write(ouf,*) 'Number of active orbitals in input  :', nacorb
      !    stop
      ! end if

      allocate(scr_cmo(ncmot),stat=ierr)
      if (ierr .ne. 0) stop 'Error allocating scr_cmo'
      !
      allocate(scr(nnashx),stat=ierr)
      if (ierr .ne. 0) stop 'Error allocating scr.1'
      !
      read(sirilu) scr_cmo  ! cmo
      read(sirilu)          ! CI coefficients
      read(sirilu) scr      ! Active 1e-density

      ! -----------
      !  CMOs
      ! -----------
      allocate(CMO(norb,norb), stat=ierr)
      if (ierr .ne. 0) stop 'read_sirifc: Error in allocation of CMO'
      CMO = 0.d0

      ij = 0
      do i = 1,norb
         do j = 1,i
            ij = ij + 1
            CMO(i,j) = scr(ij)    
            CMO(j,i) = scr(ij)    
         end do
      end do

      ! -----------
      !  1-RDM
      ! -----------
      allocate(D1MO(norb,norb), stat=ierr)
      if (ierr .ne. 0) stop 'read_sirifc: Error in allocation of D1MO'
      D1MO = 0.d0

      ! Inactive block
      do i = 1,nisht
         D1MO(i,i) = 2.d0
      end do

      ! Active block
      ij = 0
      do i = 1,nasht
         do j = 1,i
            ij = ij + 1
            ni = nisht + i
            nj = nisht + j
            D1MO(ni,nj) = scr(ij)    
            D1MO(nj,ni) = scr(ij)    
         end do
      end do
      !
      ! -----------
      !  2-RDM
      ! -----------
      deallocate(scr)
      nnashy = nnashx*nnashx
      allocate(scr(nnashy),stat=ierr)
      if (ierr .ne. 0) stop 'Error allocating scr.2'
      !
      read(sirilu)         ! Fock matrix
      read(sirilu) scr     ! Active-active 2e-density

      allocate(D2MO(norb,norb,norb,norb), stat=ierr)
      if (ierr .ne. 0) stop 'read_sirifc: Error in allocation of D2MO'

      !
      !     Inactive-Inactive
      !
      do k = 1,nisht
         do i = 1,nisht
            D2MO(i,i,k,k) = 2.d0
         end do
      end do
      !
      do i = 1,nisht
         do j = 1,nisht
            D2MO(i,j,j,i) = D2MO(i,j,j,i) - 1.d0
         end do
      end do
      !
      !     Active-Inactive
      !
      do i = 1,nisht
         do v = 1,nasht
            nv = nisht + v
            do u = 1,nasht
               nu = nisht + u
               D2MO(nu,nv,i,i) = 2.d0 * D1MO(nu,nv)
            end do
         end do
      end do
      !
      do v = 1,nasht
         nv = nisht + v
         do i = 1,nisht
            do u = 1,nasht
               nu = nisht + u
               D2MO(nu,i,i,nv) = -D1MO(nu,nv)
            end do
         end do
      end do
      !
      !     Active-Active
      !
      do y = 1,nasht
         ny = nisht + y
         do x = 1,nasht
            nx = nisht + x
            xy = indij(x,y)
            do v = 1,nasht
               nv = nisht + v
               do u = 1,nasht
                  nu = nisht + u
                  uv = indij(u,v)
                  koff = nnashx*(xy-1) + uv
                  D2MO(nu,nv,nx,ny) = 0.5d0 * scr(koff)
               end do
            end do
         end do
      end do

      !
      ! Note: do not deallocate scr
      !

   end subroutine read_sirifc

   subroutine wrimat(xmat,nrow,ncol)
      !    
      !        Subroutine to write a matrix in blocks of 5 columns 
      !    
      implicit none
      integer(kind=8) :: i, j
      integer(kind=8) :: nrow, ncol, istr, iend
      real(kind=8) :: xmat(nrow,ncol)
      integer, parameter :: maxcol = 8
      !        
      iend = 0
      do
         istr = iend + 1
         iend = iend + maxcol ; if (iend > ncol) iend = ncol
         do i = 1, nrow
            write(ouf,'(8f10.6)') (xmat(i,j),j=istr,iend)
         end do
         write(ouf,'(//)')
         if (iend == ncol) exit
      end do
      return
   end subroutine

   subroutine print_vector(label, vec, colrow, lu, fmt)
      implicit none
      character(len=*), intent(in) :: label, colrow, fmt 
      real(kind=8), dimension(:), intent(in) :: vec
      integer, intent(in) :: lu 
      !
      integer :: i

      write(lu,'(a)') label
      if (colrow.eq.'col' .or. colrow.eq.'COL') then
         ! Write it as column vector
         do i=1, size(vec,1)
            write(lu,trim(fmt)) vec(i)
         end do
      elseif (colrow.eq.'row' .or. colrow.eq.'ROW') then
         ! Write it in a row
         write(lu,trim(fmt)) ( vec(i), i=1, size(vec,1) )
      end if
      !
      return
   end subroutine print_vector 

   subroutine print_matrix(label, A, lu, fmt)
      implicit none
      character(len=*), intent(in) :: label, fmt 
      real(kind=8), dimension(:,:), intent(in) :: A
      integer, intent(in) :: lu 

      write(lu,'(a)') label
      call write_mat(A, lu, fmt, size(A, dim=2))
      write(lu,*)
      !
      return
   end subroutine print_matrix 
      
   subroutine print_energy
      E_ina = Eoe_ina + Eee_ina
      E_cross = Eoe_act + Eee_cross
      E_act = E_cross + Eee_act
      E_elec = Eoe + Eee
      E_MCSCF = E_elec + repnuc

      if (tstdbg) then
         write(ouf,'(a,2f16.10)') 'Eoe_ina, error   :', Eoe_ina, &
         & dabs(Eoe_ina-Eoe_ina_exact)
         write(ouf,'(a,2f16.10)') 'Eoe_act, error   :', Eoe_act, &
         & dabs(Eoe_act-Eoe_act_exact)
         write(ouf,'(a,2f16.10)') 'Eee_ina, error   :', Eee_ina, &
         & dabs(Eee_ina-Eee_ina_exact)
         write(ouf,'(a,2f16.10)') 'Eee_act, error   :', Eee_act, &
         & dabs(Eee_act-Eee_act_exact)
         write(ouf,'(a,2f16.10)') 'Eee_cross, error :', Eee_cross, &
         & dabs(Eee_cross-Eee_cross_exact)
         write(ouf,'(a,2f16.10)') 'Eoe, error :', Eoe, &
         & dabs(Eoe-Eoe_exact)
         write(ouf,'(a,2f16.10)') 'Eee, error :', Eee, &
         & dabs(Eee-Eee_exact)
         write(ouf,'("------")')
         write(ouf,'(a,2f16.10)') 'Inactive energy, error      :', E_ina, &
         & dabs(E_ina-siri_einactiv)
         E_cross_exact = Eoe_act_exact + Eee_cross_exact
         write(ouf,'(a,2f16.10)') 'Act-ina contribution, error :', &
         & E_cross, dabs(E_cross-E_cross_exact)
         write(ouf,'(a,2f16.10)') 'Act-act contribution, error :', Eee_act, &
         dabs(Eee_act-Eee_act_exact)
         write(ouf,'(a,2f16.10)') 'Active energy, error        :', E_act, &
         & dabs(siri_eactiv-E_act)
         write(ouf,'(a,2f16.10,/)') 'MCSCF energy, error         :', E_MCSCF, &
         & dabs(siri_emcscf-E_MCSCF)
      else
         write(ouf,'(a,2f16.10)') 'Inactive energy, error :', E_ina, &
         & dabs(E_ina-siri_einactiv)
         write(ouf,'(a,2f16.10)') 'act-ina contribution   :', E_cross
         write(ouf,'(a,2f16.10)') 'act-act contribution   :', Eee_act 
         write(ouf,'(a,2f16.10)') 'Active energy, error   :', E_act, &
         & dabs(siri_eactiv-E_act)
         write(ouf,'(a,2f16.10,/)') 'MCSCF energy, error    :', E_MCSCF, &
         & dabs(siri_emcscf-E_MCSCF)
      end if
      call flush(ouf)
   end subroutine print_energy

   subroutine output_energy(lu, label, E_MCSCF, err)
      implicit none
      integer, intent(in) :: lu 
      character(len=*), intent(in) :: label
      real(kind=8), intent(in) :: E_MCSCF, err
      !
      write(lu,'(a9,2f16.10)') label, E_MCSCF, err
      !
      return
   end subroutine output_energy 

   subroutine print_MCSCF_energy
      E_elec = Eoe + Eee
      E_MCSCF = E_elec + repnuc
      write(ouf,'(a,2f16.10)') 'MCSCF energy, error    :', E_MCSCF, &
      & dabs(siri_emcscf-E_MCSCF)
      call flush(ouf)
   end subroutine print_MCSCF_energy

end module IO
