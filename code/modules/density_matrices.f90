! +--------------------------------------------+
! | Author: Jose Antonio Quinonero Gris        |
! | Creation date: Tuesday 17:30:52 30/04/2024 |
! +--------------------------------------------+
module density_matrices
    ! 
    ! Module to contain density matrices procedures
    !
    use declarations
    use strings
    use IO
    use matmod
    !
    implicit none
    !
    contains

    subroutine calc_casscf_dalton
        write(ouf,'(a)') 'Energy calculated with Dalton algorithm for Fock matrix:'
        !
        !     1) Inactive Fock matrix
        !        iF(m,n) = h(m,n) + 2 (mn|ii) - (mi|in)
        !
        allocate(fockin(norb,norb),stat=ierr)
        if (ierr .ne. 0) stop 'calc_casscf_dalton: Error in allocation of fockin'
        !
        do m = 1,norbt
            do n = 1,norbt
                fockin(m,n) = H1MO(m,n)
                do i = 1,nisht
                    fockin(m,n) = fockin(m,n) &
                    &                     + 2.d0*H2MO(m,n,i,i) - H2MO(m,i,i,n)
                end do
            end do
        end do
        !
        !     2) Ein = h(i,i) + iF(i,i) 
        !
        eina = 0.d0
        do i = 1,nisht
            eina = eina + H1MO(i,i) + fockin(i,i)
        end do
        errina = abs(eina-siri_einactiv)
        write(ouf,'(a,2f16.10)') 'Inactive energy, error :',eina,errina
        !
        !     3) Cross contribution: Eac1 = D(u,v) iF(u,v)
        !
        eact1 = 0.d0
        do u = 1,nasht
            nu = nisht + u
            do v = 1,nasht
                nv = nisht + v
                eact1 = eact1 + D1MO(nu,nv) * fockin(nu,nv)
            end do
        end do
        !
        !     4) Q contribution: 0.5 Q(u,u)
        !        Q(p,u) = PV(uv,xy) (pv|xy)  (construct only diagonal)
        !
        allocate(qmat(nasht),stat=ierr)
        if (ierr .ne. 0) stop 'Error in allocation of qmat'
        qmat = 0.d0
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
                        !
                        koff = nnashx*(xy-1) + uv
                        qmat(u) = qmat(u) + scr(koff)*H2MO(nu,nv,nx,ny)
                    end do
                end do
            end do
        end do
        deallocate(scr)
        !
        eact2 = 0.d0
        do u = 1,nasht
            eact2 = eact2 + qmat(u)
        end do
        eact2 = 0.5d0 * eact2
        !
        eact   = eact1 + eact2
        erract = abs(siri_eactiv-eact)
        write(ouf,'(a,2f16.10)') 'Cross contribution     :',eact1
        write(ouf,'(a,2f16.10)') 'Qmat contribution      :',eact2
        write(ouf,'(a,2f16.10,/)') 'Active energy, error   :',eact,erract
        call flush(ouf)
    end subroutine calc_casscf_dalton

    subroutine calc_Eoe(D1, H1, nisht, nasht, eina, eact, Eoe)
        implicit none
        real(kind=8), dimension(:,:), intent(in) :: D1, H1 
        integer(kind=8), intent(in) :: nisht, nasht 
        real(kind=8), intent(out) :: eina, eact, Eoe
        !
        integer(kind=8) :: i, j, u, nu, v, nv 
        !
        ! Inactive energy
        !
        eina = 0.d0
        do j = 1,nisht
            do i = 1,nisht
                eina = eina + D1(i,j) * H1(i,j)
            end do
        end do
        !
        ! Active energy
        !
        eact = 0.d0
        do v = 1,nasht              ! Actually, this should go to 4)
            nv = nisht + v           ! But here it is easier to debug
            do u = 1,nasht
                nu = nisht + u
                eact = eact + D1(nu,nv) * H1(nu,nv)
            end do
        end do

        Eoe = eina + eact

        !
    end subroutine calc_Eoe

    subroutine calc_Eee(D2, H2, nisht, nasht, eina, eact, ecross, Eee)
        implicit none
        real(kind=8), dimension(:,:,:,:), intent(in) :: D2, H2 
        integer(kind=8), intent(in) :: nisht, nasht 
        real(kind=8), intent(out) :: eina, eact, ecross, Eee
        !
        integer(kind=8) :: i, j, k, l, u, nu, v, nv, x, nx, y, ny 

        ! call print_matrix('D2', D2(1,1,1:4,1:4), ouf, dfmt)

        ! Inactive energy
        !
        eina = 0.d0
        do l = 1,nisht
            do k = 1,nisht
                do j = 1,nisht
                    do i = 1,nisht
                        eina = eina + D2(i,j,k,l) * H2(i,j,k,l)
                    end do
                end do
            end do
        end do
        !
        ! Inactive-active energy
        !
        ecross = 0.d0
        do j = 1,nisht
            do i = 1,nisht
                do v = 1,nasht
                    nv = nisht + v
                    do u = 1,nasht
                        nu = nisht + u
                        ecross = ecross + D2MO(nu,nv,i,j) * H2MO(nu,nv,i,j)
                    end do
                end do
            end do
        end do
        !
        do v = 1,nasht
            nv = nisht + v
            do i = 1,nisht
                do j = 1,nisht
                    do u = 1,nasht
                        nu = nisht + u
                        ecross = ecross + D2MO(nu,j,i,nv) * H2MO(nu,j,i,nv)
                    end do
                end do
            end do
        end do
        !
        ! Active energy
        !
        eact = 0.d0
        do y = 1,nasht
            ny = nisht + y
            do x = 1,nasht
                nx = nisht + x
                do v = 1,nasht
                    nv = nisht + v
                    do u = 1,nasht
                        nu = nisht + u
                        eact = eact + D2MO(nu,nv,nx,ny)*H2MO(nu,nv,nx,ny)
                    end do
                end do
            end do
        end do

        Eee = eina + ecross + eact

        !
    end subroutine calc_Eee

    subroutine calc_Eee_sym(D2, H2, nisht, nasht, eina, eact, ecross, Eee)
        !
        ! Symmetric algorithm
        !
        implicit none
        real(kind=8), dimension(:,:,:,:), intent(in) :: D2, H2 
        integer(kind=8), intent(in) :: nisht, nasht 
        real(kind=8), intent(out) :: eina, eact, ecross, Eee
        !
        integer(kind=8) :: i, j, k, l, u, nu, v, nv, x, nx, y, ny 

        ! Inactive energy
        !
        eina = 0.d0
        do l = 1,nisht
            do k = 1,nisht
                do j = 1,nisht
                    do i = 1,nisht
                        eina = eina + D2(i,j,k,l) * H2(i,j,k,l)
                    end do
                end do
            end do
        end do
        !
        ! Inactive-active energy
        !
        ecross = 0.d0
        do j = 1,nisht
            do i = 1,nisht
                do v = 1,nasht
                    nv = nisht + v
                    do u = 1,nasht
                        nu = nisht + u
                        ecross = ecross + &
                        & D2MO(nu,nv,i,j) * H2MO(nu,nv,i,j) + &
                        & D2MO(i,j,nu,nv) * H2MO(i,j,nu,nv)
                    end do
                end do
            end do
        end do
        !
        do v = 1,nasht
            nv = nisht + v
            do i = 1,nisht
                do j = 1,nisht
                    do u = 1,nasht
                        nu = nisht + u
                        ecross = ecross + &
                        & D2MO(nu,j,i,nv) * H2MO(nu,j,i,nv) + &
                        & D2MO(j,nu,i,nv) * H2MO(j,nu,i,nv) + &
                        & D2MO(nu,j,nv,i) * H2MO(nu,j,nv,i) + &
                        & D2MO(j,nu,nv,i) * H2MO(j,nu,nv,i)
                    end do
                end do
            end do
        end do
        !
        ! Active energy
        !
        eact = 0.d0
        do y = 1,nasht
            ny = nisht + y
            do x = 1,nasht
                nx = nisht + x
                do v = 1,nasht
                    nv = nisht + v
                    do u = 1,nasht
                        nu = nisht + u
                        eact = eact + D2MO(nu,nv,nx,ny)*H2MO(nu,nv,nx,ny)
                    end do
                end do
            end do
        end do

        Eee = eina + ecross + eact

        !
    end subroutine calc_Eee_sym

    subroutine calc_Eee_SO(D2, H2, nisht, nasht, eina, eact, ecross, Eee)
        implicit none
        real(kind=8), dimension(:,:,:,:), intent(in) :: D2, H2 
        integer(kind=8), intent(in) :: nisht, nasht 
        real(kind=8), intent(out) :: eina, eact, ecross, Eee
        !
        integer(kind=8) :: nnisht, nnasht 
        integer(kind=8) :: i, j, k, l, u, nu, v, nv, x, nx, y, ny 

        ! call print_matrix('D2', D2(1,1,1:4,1:4), ouf, dfmt)

        ! For SO basis, the MO orbitals are splitted into alpha, beta. Therefore
        ! nisht <- 2*nisht ; nasht <- 2*nasht
        !
        nnisht = 2 * nisht
        nnasht = 2 * nasht

        ! Inactive energy
        !
        eina = 0.d0
        do l = 1,nnisht
            do k = 1,nnisht
                do j = 1,nnisht
                    do i = 1,nnisht
                        eina = eina + D2(i,j,k,l) * H2(i,j,k,l)
                    end do
                end do
            end do
        end do
        !
        ! Inactive-active energy
        !
        ecross = 0.d0
        do j = 1,nnisht
            do i = 1,nnisht
                do v = 1,nnasht
                    nv = nnisht + v
                    do u = 1,nnasht
                        nu = nnisht + u
                        ecross = ecross + D2MO(nu,nv,i,j) * H2MO(nu,nv,i,j)
                    end do
                end do
            end do
        end do
        !
        do v = 1,nnasht
            nv = nnisht + v
            do i = 1,nnisht
                do j = 1,nnisht
                    do u = 1,nnasht
                        nu = nnisht + u
                        ecross = ecross + D2MO(nu,j,i,nv) * H2MO(nu,j,i,nv)
                    end do
                end do
            end do
        end do
        !
        ! Active energy
        !
        eact = 0.d0
        do y = 1,nnasht
            ny = nnisht + y
            do x = 1,nnasht
                nx = nnisht + x
                do v = 1,nnasht
                    nv = nnisht + v
                    do u = 1,nnasht
                        nu = nnisht + u
                        eact = eact + D2MO(nu,nv,nx,ny)*H2MO(nu,nv,nx,ny)
                    end do
                end do
            end do
        end do

        Eee = eina + ecross + eact

        !
    end subroutine calc_Eee_SO

    subroutine calc_casscf_nonsymmetric
        write(ouf,'(a)') 'Densities algorithm (Alfredo):'
        !
        ! Inactive energy
        !
        eina = 0.d0
        do j = 1,nisht
            do i = 1,nisht
                eina = eina + D1MO(i,j) * H1MO(i,j)
            end do
        end do
        !
        do l = 1,nisht
            do k = 1,nisht
                do j = 1,nisht
                    do i = 1,nisht
                        eina = eina + D2MO(i,j,k,l) * H2MO(i,j,k,l)
                    end do
                end do
            end do
        end do
        !
        errina = abs(eina-siri_einactiv)
        write(ouf,'(a,2f16.10)') 'Inactive energy, error :',eina,errina
        call flush(ouf)
        !
        !     3) Active-inactive contribution
        !
        eact1 = 0.d0
        do v = 1,nasht              ! Actually, this should go to 4)
            nv = nisht + v           ! But here it is easier to debug
            do u = 1,nasht
                nu = nisht + u
                eact1 = eact1 + D1MO(nu,nv) * H1MO(nu,nv)
            end do
        end do
        !
        do j = 1,nisht
            do i = 1,nisht
                do v = 1,nasht
                    nv = nisht + v
                    do u = 1,nasht
                        nu = nisht + u
                        eact1 = eact1 + D2MO(nu,nv,i,j) * H2MO(nu,nv,i,j)
                    end do
                end do
            end do
        end do
        !
        do v = 1,nasht
            nv = nisht + v
            do i = 1,nisht
                do j = 1,nisht
                    do u = 1,nasht
                        nu = nisht + u
                        eact1 = eact1 + D2MO(nu,j,i,nv) * H2MO(nu,j,i,nv)
                    end do
                end do
            end do
        end do
        !
        write(ouf,'(a,2f16.10)') 'act-ina contribution   :',eact1
        call flush(ouf)
        !
        !     4) Active-active contribution
        !
        eact2 = 0.d0
        do y = 1,nasht
            ny = nisht + y
            do x = 1,nasht
                nx = nisht + x
                do v = 1,nasht
                    nv = nisht + v
                    do u = 1,nasht
                        nu = nisht + u
                        eact2 = eact2 + D2MO(nu,nv,nx,ny)*H2MO(nu,nv,nx,ny)
                    end do
                end do
            end do
        end do
        !
        write(ouf,'(a,2f16.10)') 'act-act contribution   :',eact2
        call flush(ouf)
        !
        eact   = eact1 + eact2
        erract = abs(siri_eactiv-eact)
        write(ouf,'(a,2f16.10,/)') 'Active energy, error   :',eact,erract
        call flush(ouf)
    end subroutine calc_casscf_nonsymmetric

    subroutine checknorm_D1(D1, Nelec, thres, lu)
        implicit none
        real(kind=8), dimension(:,:), intent(in) :: D1
        integer(kind=8), intent(in) :: Nelec
        real(kind=8), intent(in) :: thres 
        integer, intent(in) :: lu
        !
        if (trace(D1).lt.dble(Nelec)-thres .or. &
            & trace(D1).gt.dble(Nelec)+thres) then 
            !
            write(lu,'(a)') 'WARNING: checknorm_D1: D1 is not normalized'
            write(lu,'(a,f8.4,a,i0,/)') 'Tr[D1] = ',trace(D1),', Nelec = ',Nelec
            stop
        else
            write(lu,'(a)') 'checknorm_D1: D1 is normalized'
            write(lu,'(a,f8.4,a,i0,/)') 'Tr[D1] = ',trace(D1),', Nelec = ',Nelec
        end if
        !
        return
    end subroutine checknorm_D1

    subroutine checknorm_D2(D2, Nelec, thres, lu)
        implicit none
        real(kind=8), dimension(:,:,:,:), intent(in) :: D2
        integer(kind=8), intent(in) :: Nelec
        real(kind=8), intent(in) :: thres 
        integer, intent(in) :: lu
        !
        integer :: i, j, n 
        real(kind=8) :: suma, norma

        n = size(D2,1)

        suma = 0.d0
        do i=1, n
            do j = 1, n
                ! suma = suma + D2(i,j,i,j)  ! if (H2 in Dirac notation)
                suma = suma + D2(i,i,j,j)  ! if (H2 in Mulliken notation)
            end do
        end do

        norma = dble(nelec * (nelec - 1)) / 2.d0

        if (suma.lt.norma-thres .or. &
            & suma.gt.norma+thres) then 
            !
            write(lu,'(a)') 'WARNING: checknorm_D2: D2 is not normalized'
            write(lu,'(a,f8.4,a,f8.4,/)') 'sum_{pq} D^(2)_{pqpq} = ',suma,', N(N-1)/2 = ',norma
            ! TODO: descomentar el stop
            ! stop
        else
            write(lu,'(a)') 'checknorm_D2: D2 is normalized'
            write(lu,'(a,f8.4,a,f8.4,/)') 'sum_{pq} D^(2)_{pqpq} = ',suma,', N(N-1)/2 = ',norma
        end if
        !
        return
    end subroutine checknorm_D2

    subroutine symmetrize_D2(D2, nisht, nocct, tstdbg, lu)
        !
        ! Subroutine to symmetrize the 2-RDM
        !
        implicit none
        real(kind=8), dimension(:,:,:,:), intent(inout) :: D2
        integer(kind=8), intent(in) :: nisht, nocct
        logical, intent(in) :: tstdbg
        integer, intent(in)  :: lu 
        !
        logical :: sym1, sym12, sym2, symtot
        integer(kind=8) :: i, j, k, l
        real(kind=8) :: x1, x2 

        ! D2(j,i,l,k) = D2(i,j,k,l)   ! Hermiticity
        ! D2(k,l,i,j) = D2(i,j,k,l)   ! Pair exchange
        ! D2(j,i,k,l) = -D2(i,j,k,l)  ! Antisymmetry for change
        ! D2(i,j,l,k) = -D2(i,j,k,l)  ! of intern indexes

        !
        !     a. Inactive-Inactive
        !
        do l = 1,nisht
            do k = 1,nisht
                do j = 1,nisht
                    do i = 1,j-1
                        x1 = D2(i,j,k,l)
                        x2 = D2(j,i,k,l)
                        D2(i,j,k,l) = 0.5d0 * (x1+x2)
                        D2(j,i,k,l) = 0.5d0 * (x1+x2)
                    end do
                end do
            end do
        end do
        !
        do l = 1,nisht
            do k = 1,l-1
                do j = 1,nisht
                    do i = 1,nisht
                        x1 = D2(i,j,k,l)
                        x2 = D2(i,j,l,k)
                        D2(i,j,k,l) = 0.5d0 * (x1+x2)
                        D2(i,j,l,k) = 0.5d0 * (x1+x2)
                    end do
                end do
            end do
        end do

        sym1  = .true.
        sym2  = .true.
        sym12 = .true.
        !
        do l = 1,nisht
            do k = 1,nisht
                do j = 1,nisht
                    do i = 1,nisht
                        if (D2(i,j,k,l) .ne. D2(j,i,k,l)) sym1  = .false.
                        if (D2(i,j,k,l) .ne. D2(i,j,l,k)) sym2  = .false.
                        if (D2(i,j,k,l) .ne. D2(k,l,i,j)) sym12 = .false.
                        ! write(6,'(4i3,4f16.10)') i,j,k,l,D2(i,j,k,l),
                        ! & D2(j,i,k,l),D2(i,j,l,k),D2(k,l,i,j)
                    end do
                end do
            end do
        end do
        !
        if (tstdbg) then 
            write(lu,'(/,a)') 'Symmetry of inactive-inactive block of D2'
            write(lu,'(a,l2)') ' i-j  swapping :',sym1
            write(lu,'(a,l2)') ' k-l  swapping :',sym2
            write(lu,'(a,l2)') 'ij-kl swapping :',sym12
        end if
        !
        symtot = sym1 .and. sym2 .and. sym12
        if (.not. symtot) stop 'symmetrize_D2: Gamma_ijkl non-symmetric'
        !
        !     b. Active-Inactive Coulomb (uv,ij)
        !
        do j = 1,nisht
            do i = 1,nisht
                do v = nisht+1,nocct
                    do u = nisht+1,nocct
                        x1 = D2(u,v,i,j) 
                        x2 = D2(i,j,u,v)
                        D2(u,v,i,j) = 0.5d0 * (x1+x2)
                        D2(i,j,u,v) = 0.5d0 * (x1+x2)
                    end do
                end do
            end do
        end do
        !
        sym1  = .true.
        sym2  = .true.
        sym12 = .true.
        !
        do j = 1,nisht
            do i = 1,nisht
                do v = nisht+1,nocct
                    do u = nisht+1,nocct
                        if (D2(u,v,i,j) .ne. D2(v,u,i,j)) sym1  = .false.
                        if (D2(u,v,i,j) .ne. D2(u,v,j,i)) sym2  = .false.
                        if (D2(u,v,i,j) .ne. D2(i,j,u,v)) sym12 = .false.
                    end do
                end do
            end do
        end do
        !
        if (tstdbg) then
            write(lu,'(/,a)') 'Symmetry of act-inact Coulomb block of D2'
            write(lu,'(a,l2)') ' u-v  swapping :',sym1
            write(lu,'(a,l2)') ' i-j  swapping :',sym2
            write(lu,'(a,l2)') 'uv-ij swapping :',sym12
        end if
        !
        symtot = sym1 .and. sym2 .and. sym12
        if (.not. symtot) stop 'symmetrize_D2: Gamma_uvij non-symmetric'
        !
        !     c. Active-Inactive Exchange (uj,iv)
        !
        do v = nisht+1, nocct
            do i = 1,nisht
                do j = 1,nisht
                    do u = nisht+1,nocct
                        x1 = D2(u,j,i,v)
                        x2 = D2(j,u,i,v)
                        D2(u,j,i,v) = 0.5d0 * (x1+x2)
                        D2(j,u,i,v) = 0.5d0 * (x1+x2)
                    end do
                end do
            end do
        end do
        !
        do v = nisht+1, nocct
            do i = 1,nisht
                do j = 1,nisht
                    do u = nisht+1,nocct
                        x1 = D2(u,j,i,v)
                        x2 = D2(u,j,v,i)
                        D2(u,j,i,v) = 0.5d0 * (x1+x2)
                        D2(u,j,v,i) = 0.5d0 * (x1+x2)
                    end do
                end do
            end do
        end do
        !
        sym1  = .true.
        sym2  = .true.
        sym12 = .true.
        !
        do j = 1,nisht
            do i = 1,nisht
                do v = nisht+1,nocct
                    do u = nisht+1,nocct
                        if (D2(u,v,i,j) .ne. D2(v,u,i,j)) sym1  = .false.
                        if (D2(u,v,i,j) .ne. D2(u,v,j,i)) sym2  = .false.
                        if (D2(u,v,i,j) .ne. D2(i,j,u,v)) sym12 = .false.
                    end do
                end do
            end do
        end do
        !
        if (tstdbg) then
            write(lu,'(/,a)') 'Symmetry of act-inact Exchange block of D2'
            write(lu,'(a,l2)') ' u-v  swapping :',sym1
            write(lu,'(a,l2)') ' i-j  swapping :',sym2
            write(lu,'(a,l2)') 'uv-ij swapping :',sym12
        end if
        !
        symtot = sym1 .and. sym2 .and. sym12
        if (.not. symtot) stop 'symmetrize_D2: Gamma_ujiv non-symmetric'
        !
        !     d. Active-Active
        !
        sym1  = .true.
        sym2  = .true.
        sym12 = .true.
        !
        do y = nisht+1,nocct
            do x = nisht+1,nocct
                do v = nisht+1,nocct
                    do u = nisht+1,nocct
                        if (D2(u,v,x,y) .ne. D2(v,u,x,y)) sym1  = .false.
                        if (D2(u,v,x,y) .ne. D2(u,v,y,x)) sym2  = .false.
                        if (D2(u,v,x,y) .ne. D2(x,y,u,v)) sym12 = .false.
                    end do
                end do
            end do
        end do
        !
        if (tstdbg) then
            write(lu,'(/,a)') 'Symmetry of active-active block of D2'
            write(lu,'(a,l2)') ' u-v  swapping :',sym1
            write(lu,'(a,l2)') ' x-y  swapping :',sym2
            write(lu,'(a,l2)') 'uv-xy swapping :',sym12
        end if
        !
        symtot = sym1 .and. sym2 .and. sym12
        if (.not. symtot) stop 'symmetrize_D2: Gamma_uvxy non-symmetric'

        !
        return
    end subroutine symmetrize_D2 

    subroutine check_virtorbs_D1(D1, nocct, norbt, ouf)
        !
        ! Subroutine to check that all virtual orbitals of D^(1) are zero
        !
        implicit none
        real(kind=8), dimension(:,:), intent(in) :: D1
        integer(kind=8), intent(in) :: nocct, norbt
        integer, intent(in) :: ouf 
        !
        integer(kind=8) :: i, j

        ! Check D^(1)
        do j = nocct+1,norbt
            do i = nocct+1,norbt
                if (D1(i,j) .ne. 0.d0) then
                    write(ouf,'(a,1x,i0,",",i0,1x,a)') &
                    & 'WARNING check_virtorbs: D^(1) with virtual index', i, j,&
                    'non-zero'
                    stop
                end if
            end do
        end do
        !
        return
    end subroutine check_virtorbs_D1

    subroutine check_virtorbs_D2(D2, nocct, norbt, ouf)
        !
        ! Subroutine to check that all virtual orbitals of D^(2) are zero
        !
        implicit none
        real(kind=8), dimension(:,:,:,:), intent(in) :: D2
        integer(kind=8), intent(in) :: nocct, norbt
        integer, intent(in) :: ouf 
        !
        integer(kind=8) :: i, j, k, l

        ! Check D^(2)
        do l = nocct+1,norbt
            do k = nocct+1,norbt
                do j = nocct+1,norbt
                    do i = nocct+1,norbt
                        if (D2(i,j,k,l) .ne. 0.d0) then
                            write(ouf,'(a,1x,i0,",",i0,1x,a)') &
                            &'WARNING check_virtorbs: D^(2) with virtual index',&
                            & i, j, 'non-zero'
                            stop
                        end if
                    end do
                end do
            end do
        end do
        
        !
        return
    end subroutine check_virtorbs_D2

    subroutine check_SO_ON(ON)
        !
        ! Subroutine to check that ON is in spin-orbital basis, i.e.
        ! 0 <= n_p <= 1 forall p
        !
        implicit none
        real(kind=8), dimension(:), intent(in) :: ON
        integer :: i 

        do i = 1, size(ON)
            if (ON(i).gt.1.d0) then
                stop 'check_SO_ON: Error: ON(i) > 1'
            end if
        end do

        return
    end subroutine check_SO_ON 

    subroutine MO_to_SO_D1(AMO, ASO)
        !
        ! Subroutine to transform 1-RDM from MO basis to spin-orbital basis, i.e.
        ! 0 <= nu_p <= 2 forall p in MOs -> 0 <= n_p <= 1 forall p in SOs
        !
        implicit none
        real(kind=8), dimension(:,:), intent(in) :: AMO
        real(kind=8), dimension(:,:), allocatable, intent(out) :: ASO
        !
        integer :: ierr, i, j, n

        ! Get size of the matrix & asure it is square
        n = size(AMO, dim=1)
        if (n .ne. size(AMO, dim=2)) &
            & stop 'MO_to_SO_D1: Error: input matrix not square'

        ! Allocate ASO
        allocate(ASO(2*n, 2*n), stat=ierr)
        if (ierr .ne. 0) stop 'MO_to_SO_D1: Error in allocation of ASO'
        ASO = 0.d0

        ! Transform to spin-orbital basis as:
        ! ASO(ialpha, jalpha) = 1/2 AMO(i, j)
        ! ASO(ibeta, jbeta) = 1/2 AMO(i, j)
        ! ASO(ialpha, jbeta) = ASO(ibeta, jalpha) = 0.d0
        do i = 1, n
            do j = 1, n
                ASO(2*i-1, 2*j-1) = 0.5d0 * AMO(i, j)  ! alpha-alpha
                ASO(2*i, 2*j) = 0.5d0 * AMO(i, j)      ! beta-beta
            end do
        end do

        !
        return
    end subroutine MO_to_SO_D1 

    subroutine MO_to_SO_D2(AMO, ASO)
        !
        ! Subroutine to transform 2-RDM from MO basis to spin-orbital basis, i.e.
        ! 0 <= nu_p <= 2 forall p in MOs -> 0 <= n_p <= 1 forall p in SOs
        !
        implicit none
        real(kind=8), dimension(:,:,:,:), intent(in) :: AMO
        real(kind=8), dimension(:,:,:,:), allocatable, intent(out) :: ASO
        !
        integer :: ierr, i, j, k, l, n

        ! Get size of the matrix & asure it is square
        n = size(AMO, dim=1)
        if (n .ne. size(AMO, dim=2) .or. &
            & n .ne. size(AMO, dim=3) .or. &
            & n .ne. size(AMO, dim=4) ) &
            & stop 'MO_to_SO_D2: Error: input matrix not square'

        ! Allocate ASO
        allocate(ASO(2*n, 2*n, 2*n, 2*n), stat=ierr)
        if (ierr .ne. 0) stop 'MO_to_SO_D2: Error in allocation of ASO'
        ASO = 0.d0

        ! Transform to spin-orbital basis as:
        ! The rest = 0
        do l = 1, n
            do k = 1, n
                do j = 1, n
                    do i = 1, n
                        ! alpha alpha alpha alpha
                        ASO(2*i-1, 2*j-1, 2*k-1, 2*l-1) = 0.5d0 * AMO(i,j,k,l)
                        ! beta beta beta beta
                        ASO(2*i, 2*j, 2*k, 2*l) = 0.5d0 * AMO(i,j,k,l)
                        ! alpha alpha beta beta
                        ASO(2*i-1, 2*j-1, 2*k, 2*l) = 0.5d0 * AMO(i,j,k,l)
                        ! beta beta alpha alpha
                        ASO(2*i, 2*j, 2*k-1, 2*l-1) = 0.5d0 * AMO(i,j,k,l)
                    end do
                end do
            end do
        end do

        !
        return
    end subroutine MO_to_SO_D2 

    subroutine MO_to_SO_H1(H1MO, H1SO)
        !
        ! Subroutine to transform one-electron integrals in spacial orbital
        ! basis (p|h|q) to spin-orbital basis <p|h|q>
        !
        implicit none
        real(kind=8), dimension(:,:), intent(in) :: H1MO
        real(kind=8), dimension(:,:), allocatable, intent(out) :: H1SO
        !
        integer :: p, q, ierr, n 

        ! Check for inconsistencies in the dimensions
        n = size(H1MO, dim=1)
        if (n.ne.size(H1MO, dim=2)) &
            & stop 'MO_to_SO_H1: Error: Inconsistent dimensions'

        allocate(H1SO(2*n, 2*n), stat=ierr)
        if (ierr .ne. 0) stop 'MO_to_SO_H1: Error in allocation of H1SO'
        H1SO = 0.d0

        do p = 1, n
            do q = 1, n
                ! alpha-alpha spin-orbitals
                H1SO(2*p-1, 2*q-1) = H1MO(p, q)
                ! beta-beta spin-orbitals
                H1SO(2*p, 2*q) = H1MO(p, q)
                ! mixed alpha-beta -> 0
            end do
        end do
        !
        return
    end subroutine MO_to_SO_H1 

    subroutine MO_to_SO_H2(H2MO, H2SO)
        !
        ! Subroutine to transform two-electron integrals in spacial orbital
        ! basis <pq|rs> to spin-orbital basis <pq|rs>
        !
        implicit none
        real(kind=8), dimension(:,:,:,:), intent(in) :: H2MO
        real(kind=8), dimension(:,:,:,:), allocatable, intent(out) :: H2SO
        !
        integer :: p, q, r, s, ierr, n 

        ! Check for inconsistencies in the dimensions
        n = size(H2MO, dim=1)
        if (n.ne.size(H2MO, dim=2) .or. &
            &   n.ne.size(H2MO, dim=3) .or. &
            &   n.ne.size(H2MO, dim=4)) &
            & stop 'MO_to_SO_H2: Error: Inconsistent dimensions'

        allocate(H2SO(2*n, 2*n, 2*n, 2*n), stat=ierr)
        if (ierr .ne. 0) stop 'MO_to_SO_H2: Error in allocation of H2SO'
        H2SO = 0.d0

        do p = 1, n
            do q = 1, n
                do r = 1, n
                    do s = 1, n

                        ! alpha alpha alpha alpha
                        H2SO(2*p-1, 2*q-1, 2*r-1, 2*s-1) = H2MO(p, q, r, s)
                        ! beta beta beta beta
                        H2SO(2*p, 2*q, 2*r, 2*s) = H2MO(p, q, r, s)
                        ! ! alpha alpha beta beta
                        ! H2SO(2*p-1, 2*q-1, 2*r, 2*s) = H2MO(p, q, r, s)
                        ! ! beta beta alpha alpha 
                        ! H2SO(2*p, 2*q, 2*r-1, 2*s-1) = H2MO(p, q, r, s)

                        ! mixed alpha-beta -> 0
                    end do
                end do
            end do
        end do
        !
        return
    end subroutine MO_to_SO_H2 

    subroutine H1D1_energy(D1, H1, E1)
        !
        ! Compute the non-ee part of the energy: tr[D^(1) H^(1)]
        !
        ! D1: one-electron density matrix
        ! H1: one-electron integral matrix
        !
        ! E1: energy
        !
        implicit none
        real(kind=8), dimension(:,:), intent(in) :: D1, H1
        real(kind=8), intent(out) :: E1
        !
        integer :: i, j, n 

        ! Get the dimensions
        n = size(D1,1)
        if (n.ne.size(D1,2) .or. n.ne.size(H1,1) .or. n.ne.size(H1,2)) &
            & stop 'H1D1_energy: Error: Inconsistent dimensions'

        ! Compute the sum
        E1 =  0.d0
        do j = 1, n
            do i = 1, n
                E1 = E1 + H1(i,j) * D1(i,j)
            end do
        end do

        !
        return
    end subroutine H1D1_energy 

    subroutine Eee_D2(D2, H2, E2)
        !
        ! Compute the e-e part of the energy from the exact D^(2):
        ! EeeD2 = 1/2 sum_{pqrs} (pq|rs) D2_{pqrs}
        !
        ! D2: two-electron density matrix
        ! H2: two-electron integral matrix
        !
        ! E2: energy (value of Eee)
        !
        implicit none
        real(kind=8), dimension(:,:,:,:), intent(in) :: D2, H2
        real(kind=8), intent(out) :: E2
        !
        integer :: i, j, k, l, n, n1, n2, n3, n4, m1, m2, m3, m4

        ! Get the dimensions
        n1 = size(D2,1) ; n2 = size(D2,2) ; n3 = size(D2,3) ; n4 = size(D2,4)
        m1 = size(H2,1) ; m2 = size(H2,2) ; m3 = size(H2,3) ; m4 = size(H2,4)
        if (n1.ne.m1 .or. n2.ne.m2 .or. n3.ne.m3 .or. n4.ne.m4) &
            & stop 'HF_energy: Error: Inconsistent dimensions'
        n = n1

        ! Compute the sum
        E2 =  0.d0
        do l = 1, n
            do k = 1, n
                do j = 1, n
                    do i = 1, n
                        ! TODO
                        E2 = E2 + D2(i,j,k,l) * H2(k,l,i,j)  ! se supone que es asi
                        ! E2 = E2 + H2(i,j,k,l) * D2(i,j,k,l) ! pero asi da energias mas cercanas a las aproxs
                    end do
                end do
            end do
        end do

        ! Compute the energy
        E2 = 0.5d0 * E2

        !
        return
    end subroutine Eee_D2 

    subroutine energy(E1, E2, Etot)
        !
        ! Subroutine to compute the total energy as a function of E_oe, E_ee:
        !
        implicit none
        real(kind=8), intent(in) :: E1, E2
        real(kind=8), intent(out) :: Etot
        !
        Etot = E1 + E2
        !
        return
    end subroutine energy 

    subroutine IA_energy(D1, D2, H1, H2, ninac, nac, Einac, Eac)
        !
        ! Subroutine to compute the inactive and active energies
        !
        implicit none
        real(kind=8), dimension(:,:), intent(in) :: D1, H1
        real(kind=8), dimension(:,:,:,:), intent(in) :: D2, H2
        integer, intent(in) :: ninac, nac
        real(kind=8), intent(out) :: Einac, Eac
        !
        integer :: n, ierr 
        real(kind=8) :: E1inac, E2inac, E1ac, E2ac
        real(kind=8), dimension(:,:), allocatable :: D1inac, H1inac, D1ac, H1ac 
        real(kind=8), dimension(:,:,:,:), allocatable :: D2inac, H2inac, D2ac, H2ac 

        !
        ! Extract the inactive block
        !
        allocate(D1inac(ninac,ninac), stat=ierr)
        if (ierr .ne. 0) stop 'IA_energy: Error in allocation of D1inac'
        !
        allocate(H1inac(ninac,ninac), stat=ierr)
        if (ierr .ne. 0) stop 'IA_energy: Error in allocation of H1inac'
        !
        allocate(D2inac(ninac,ninac,ninac,ninac), stat=ierr)
        if (ierr .ne. 0) stop 'IA_energy: Error in allocation of D2inac'
        !
        allocate(H2inac(ninac,ninac,ninac,ninac), stat=ierr)
        if (ierr .ne. 0) stop 'IA_energy: Error in allocation of H2inac'
        !
        D1inac = D1(1:ninac, 1:ninac)
        H1inac = H1(1:ninac, 1:ninac)
        D2inac = D2(1:ninac, 1:ninac, 1:ninac, 1:ninac)
        H2inac = H2(1:ninac, 1:ninac, 1:ninac, 1:ninac)

        !
        ! Extract the active block
        !
        allocate(D1ac(nac,nac), stat=ierr)
        if (ierr .ne. 0) stop 'IA_energy: Error in allocation of D1ac'
        !
        allocate(H1ac(nac,nac), stat=ierr)
        if (ierr .ne. 0) stop 'IA_energy: Error in allocation of H1ac'
        !
        allocate(D2ac(nac,nac,nac,nac), stat=ierr)
        if (ierr .ne. 0) stop 'IA_energy: Error in allocation of D2ac'
        !
        allocate(H2ac(nac,nac,nac,nac), stat=ierr)
        if (ierr .ne. 0) stop 'IA_energy: Error in allocation of H2ac'
        !
        n = ninac + nac
        D1ac = D1(ninac+1:n, ninac+1:n)
        H1ac = H1(ninac+1:n, ninac+1:n)
        D2ac = D2(ninac+1:n, ninac+1:n, ninac+1:n, ninac+1:n)
        H2ac = H2(ninac+1:n, ninac+1:n, ninac+1:n, ninac+1:n)

        !
        ! Inactive energy:
        !
        ! Compute E_oe
        call H1D1_energy(D1inac, H1inac, E1inac)
        ! Compute E_ee
        call Eee_D2(D2inac, H2inac, E2inac)
        ! Compute total energy
        call energy(E1inac, E2inac, Einac)

        !
        ! Active energy:
        !
        ! Compute E_oe
        call H1D1_energy(D1ac, H1ac, E1ac)
        ! Compute E_ee
        call Eee_D2(D2ac, H2ac, E2ac)
        ! Compute total energy
        call energy(E1ac, E2ac, Eac)

        !
        return
    end subroutine IA_energy

    subroutine Minkowski_distance_4D(p, D2, D2p, d)
        !
        ! Subroutine to compute the Minkowski distance between two 4D matrices,
        ! treated as high-dimensional vectors.
        ! Each element of the 4D matrix is considered as a coordinate in this
        ! 4D space.
        !
        ! The distance of order p between two points (x,y) in n-dimensional
        ! space is given by
        ! d(x,y) = (sum_{i=1}^n (x_i - y_i)^p)^(1/p)
        !
        ! - p=1   -> Manhattan distance
        ! - p=2   -> Euclidean distance
        ! - p->oo -> Chebyshev distance
        !
        implicit none
        integer, intent(in) :: p 
        real(kind=8), dimension(:,:,:,:), intent(in) :: D2, D2p
        real(kind=8), intent(out) :: d
        !
        integer :: i, j, k, l, n
        integer :: n1, n2, n3, n4, np1, np2, np3, np4 

        ! Check for inconsistencies in the dimension
        n1 = size(D2, dim=1) ; n2 = size(D2, dim=2)
        n3 = size(D2, dim=3) ; n4 = size(D2, dim=4)
        np1 = size(D2p, dim=1) ; np2 = size(D2p, dim=2)
        np3 = size(D2p, dim=3) ; np4 = size(D2p, dim=4)
        n = n1
        if (n.ne.n2 .or. n.ne.n3 .or. n.ne.n4 .or. &
            &n.ne.np1 .or. n.ne.np2 .or. n.ne.np3 .or. n.ne.np4) &
            & stop 'ERROR Minkowski_distance_4D: Inconsistent dimensions'

        ! Compute the sum
        d = 0.d0
        do l = 1, n
            do k = 1, n
                do j = 1, n
                    do i = 1, n
                        d = d + dabs(D2(i,j,k,l) - D2p(i,j,k,l))**p
                    end do
                end do
            end do
        end do

        ! Compute the distance
        d = d**(1.d0/dble(p))

        !
        return
    end subroutine Minkowski_distance_4D 

    subroutine D1_from_D2(D2, nelec, D1)
        !
        ! Reconstruct D^(1) form D^(2) as
        ! D1_{ik} = 2/(N-1) sum_j D2_{ijkj}
        !
        implicit none
        real(kind=8), dimension(:,:,:,:), intent(in) :: D2
        integer, intent(in) :: nelec 
        real(kind=8), dimension(:,:), allocatable, intent(out) :: D1
        !
        integer :: i, j, k, n, ierr
        real(kind=8) :: suma 

        n = size(D2, dim=1)
        allocate(D1(n,n), stat=ierr)
        if (ierr .ne. 0) stop 'D1_from_D2: Error in allocation of D1'

        suma = 0.d0
        do i = 1, n
            do k = 1, n
                do j = 1, n
                    suma = suma + D2(i,j,k,j)
                end do
                D1(i,k) = 2.d0/dble(nelec-1) * suma
            end do
        end do

        !
        return
    end subroutine D1_from_D2 

end module density_matrices 
