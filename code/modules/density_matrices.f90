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
    
    subroutine checknorm_D1(D1, Nelec, thres, lu)
        implicit none
        real(kind=8), dimension(:,:), intent(in) :: D1
        integer, intent(in) :: Nelec, lu
        real(kind=8), intent(in) :: thres 
        !
        if (trace(D1).lt.dble(Nelec)-thres .or. &
          & trace(D1).gt.dble(Nelec)+thres) then 
            !
            write(lu,'(a)') 'checknorm_D1: D1 is not normalized'
            write(lu,'(a,f18.4,a,i0)') 'Tr[D1] = ',trace(D1),', Nelec = ',Nelec
            stop
        end if
        !
        return
    end subroutine checknorm_D1

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

        ! TODO: confirmar que esta bien
        ! Transform to spin-orbital basis as:
        ! ASO(ialpha, jalpha, kalpha, lalpha) = 1/2 AMO(i, j, k, l)
        ! ASO(ibeta, jbeta, kbeta, lbeta) = 1/2 AMO(i, j, k, l)
        ! The rest = 0
        do l = 1, n
            do k = 1, n
                do j = 1, n
                    do i = 1, n
                        ASO(2*i-1, 2*j-1, 2*k-1, 2*l-1) = 0.5d0 * AMO(i,j,k,l)
                        ASO(2*i, 2*j, 2*k, 2*l) = 0.5d0 * AMO(i,j,k,l)
                    end do
                end do
            end do
        end do
        
        !
        return
    end subroutine MO_to_SO_D2 

    ! subroutine MO_to_SO_H2(H2MO, H2SO)
    !     implicit none
    !     real(kind=8), dimension(:,:,:,:), intent(in) :: H2MO
    !     real(kind=8), dimension(:,:,:,:), allocatable, intent(out) :: H2SO
    !     !
    !     integer :: p, q, r, s, ps, qs, rs, ss, ierr, n 
    !
    !     ! Check for inconsistencies in the dimensions
    !     n = size(H2MO, dim=1)
    !     if (n.ne.size(H2MO, dim=2) .or. &
    !     &   n.ne.size(H2MO, dim=3) .or. &
    !     &   n.ne.size(H2MO, dim=4)) &
    !     & stop 'MO_to_SO_H2: Error: Inconsistent dimensions'
    !
    !     allocate(H2SO(2*n, 2*n, 2*n, 2*n), stat=ierr)
    !     if (ierr .ne. 0) stop 'MO_to_SO_H2: Error in allocation of H2SO'
    !     H2SO = 0.d0
    !
    !     do p = 1, n
    !         do q = 1, n
    !             do r = 1, n
    !                 do s = 1, n
    !
    !                     ! alpha-alpha spin-orbitals
    !                     ps = 2*p - 1
    !                     qs = 2*q - 1
    !                     rs = 2*r - 1
    !                     ss = 2*s - 1
    !                     H2SO(ps, qs, rs, ss) = H2MO(p, q, r, s)
    !
    !                     ! beta-beta spin-orbitals
    !                     ps = 2*p
    !                     qs = 2*q
    !                     rs = 2*r
    !                     ss = 2*s
    !                     H2SO(ps, qs, rs, ss) = H2MO(p, q, r, s)
    !
    !                     ! mixed alpha-beta -> 0
    !                 end do
    !             end do
    !         end do
    !     end do
    !     !
    !     return
    ! end subroutine MO_to_SO_H2 

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
                        E2 = E2 + H2(i,j,k,l) * D2(i,j,k,l)
                    end do
                end do
            end do
        end do

        ! Compute the energy
        E2 = 0.5d0 * E2
        
        !
        return
    end subroutine Eee_D2 

    subroutine Eee_ELS(ONin, H2, Eee)
        implicit none
        real(kind=8), dimension(:), intent(in) :: ONin
        real(kind=8), dimension(:,:,:,:), intent(in) :: H2
        real(kind=8), intent(out) :: Eee
        !
        integer :: i, j, n, ierr 
        real(kind=8) :: sum1 
        real(kind=8), dimension(:), allocatable :: ON

        ! Check for inconsistencies in dimension
        n = size(ONin)
        if (n.ne.size(H2, dim=1) .or. &
          & n.ne.size(H2, dim=2) .or. &
          & n.ne.size(H2, dim=3) .or. &
          & n.ne.size(H2, dim=4)) &
          stop 'Eee_ELS: Error: inconsistent dimensions'

        ! Check that ON is in spin-orbital basis
        ! call check_SO_ON(ON)

        ! Check that ON is sorted in descending order
        allocate(ON(n), stat=ierr)
        if (ierr .ne. 0) stop 'Eee_ELS: Error in allocation of ON'
        ON = ONin
        call asure_descending(ONin, ON)

        ! Compute the sum
        sum1 = 0.d0
        do i = 1, n
            do j = 1, n
                sum1 = sum1 + GLS(i,j) * H2(i,j,j,i)
            end do
        end do

        ! Compute Eee
        Eee = 0.5d0 * sum1

        ! ---------- Intern functions ---------- 
        contains
        function GLS(p,q) result(res)
            implicit none
            integer, intent(in) :: p,q
            real(kind=8) :: res 
        
            if (p.eq.q) then 
                res = ON(p)
            elseif (p.eq.1 .and. q.gt.1 .or. p.gt.1 .and. q.eq.1) then 
                res = - dsqrt(ON(p) * ON(q))
            else
                res = dsqrt(ON(p) * ON(q))
            end if
        
        end function GLS 
    end subroutine Eee_ELS 

    subroutine Eee_BBC(level, ON, H2, Eee)
        implicit none
        character(len=*), intent(in) :: level 
        real(kind=8), dimension(:), intent(in) :: ON
        real(kind=8), dimension(:,:,:,:), intent(in) :: H2
        real(kind=8), intent(out) :: Eee
        !
        integer :: i, j, ierr, n
        real(kind=8) :: sum1, sum2, ni, nj
        real(kind=8), dimension(:,:), allocatable :: GBBC 

        ! Check for inconsistencies in the dimensions
        n = size(ON)
        if (n.ne.size(H2, dim=1) .or. &
        &   n.ne.size(H2, dim=2) .or. &
        &   n.ne.size(H2, dim=3) .or. &
        &   n.ne.size(H2, dim=4)) &
        & stop 'Eee_BBC: Error: Inconsistent dimensions'

        ! Check that ON is in spin-orbital basis
        ! call check_SO_ON(ON)

        ! Compute the GBBC matrix
        allocate(GBBC(n,n), stat=ierr)
        if (ierr .ne. 0) stop 'Eee_BBC: Error in allocation of GBBC'
        GBBC = 0.d0

        ! Compute the GBBC matrix
        call calc_GBBC

        ! Compute first sum
        sum1 = 0.d0
        do i = 1, n
            do j = 1, n
                sum1 = sum1 + ON(i) * ON(j) * H2(i,i,j,j)
            end do
        end do

        ! Compute second sum
        sum2 = 0.d0
        do i = 1, n
            do j = 1, n
                sum2 = sum2 + GBBC(i,j) * H2(i,j,j,i)
            end do
        end do

        ! Compute the energy
        Eee = 0.5d0 * sum1 + 0.5d0 * sum2
        
        ! ---------- Intern functions ---------- 
        contains
        !
        function is_weak(np) result(isweak)
            implicit none
            real(kind=8), intent(in) :: np
            logical :: isweak
            !
            if (np.lt.0.5d0) then
                isweak = .true.
            else
                isweak = .false.
            end if
        endfunction is_weak

        function is_frontier(np) result(isfrontier)
            implicit none
            real(kind=8), intent(in) :: np
            logical :: isfrontier
            ! TODO: comprobar que esta definicion de frontier orbital esta bien
            if (np.gt.0.5d0 .and. np.lt.1.0d0) then
                isfrontier = .true.
            else
                isfrontier = .false.
            end if
        endfunction is_frontier

        subroutine calc_GBBC
            if (level.eq.'1' .or. level.eq.'BBC1') then
                ! Corrected BB functional, BBC1
                do i = 1, n
                    do j = 1, n
                        ni = ON(i) ; nj = ON(j)
                        if (j.ne.i .and. is_weak(ni) .and. is_weak(nj)) then
                            ! For p/=q and p,q are weakly coupled (np,nq < 1/2)
                            GBBC(i,j) = dsqrt(ni * nj)
                        else
                            GBBC(i,j) = - dsqrt(ni * nj)
                        end if
                    end do 
                end do
            elseif (level.eq.'3' .or. level.eq.'BBC3') then
                do i = 1, n
                    do j = 1, n
                        ni = ON(i) ; nj = ON(j)
                        if ((i.ne.j .and. is_weak(ni) .and. is_weak(nj)) .or. &
                            ! For p/=q and p,q weak
                            &   (is_weak(ni) .and. is_frontier(nj)) .or. &
                            ! Or p weak, q frontier (weak) 
                            &   (is_frontier(ni) .and. is_weak(nj))) then
                            ! Or p frontier (weak), q weak 
                            GBBC(i, j) = dsqrt(ni * nj)
                        elseif ((i.ne.j .and. .not.is_weak(nj) .and. &
                            & .not.is_weak(nj)) .or. &
                            ! For p/=q and p,q strong
                            &   (.not.is_weak(ni) .and. is_frontier(nj)) .or. &
                            ! Or p strong, q frontier
                            &   (is_frontier(ni) .and. .not.is_weak(nj))) then
                            ! Or p frontier, q strong 
                            GBBC(i, j) = -ni * nj
                        elseif (i.eq.j .and. .not.is_frontier(ni)) then
                            ! For p=q and p is not frontier
                            GBBC(i, j) = -ni**2
                        else
                            !  Otherwise
                            GBBC(i, j) = -dsqrt(ni * nj)
                        end if
                    end do 
                end do
            else
                stop 'Eee_BBC: Error: Invalid level. Available: BBC1, BBC3'
            end if
        end subroutine calc_GBBC 
    end subroutine Eee_BBC 

    ! subroutine Eee_PNOF(ON, thres)
    !     implicit none
    !     real(kind=8), dimension(:), intent(in) :: ON
    !     real(kind=8), intent(in) :: thres 
    !     !
    !     integer :: ierr, i, j, n, occON
    !     real(kind=8), dimension(:), allocatable :: ONaux, ONred
    !     real(kind=8), dimension(:,:), allocatable :: pairs
    !
    !     !
    !     ! For PNOF5:
    !     !
    !     ! Extract the density matrix only with occupied orbitals, i.e.
    !     ! D'_ij <- D^(1)_ij > 0
    !     ! which, for the ON vector read
    !     ! ON'_i <- ON_i > 0
    !     ! Get dimensions of D1
    !     n = size(ON)
    !     allocate(ONaux(n), stat=ierr)
    !     if (ierr .ne. 0) stop 'Eee_PNOF: Error in allocation of ONaux'
    !     ONaux = 0.d0
    !     occON = 1
    !     do i = 1, n
    !         if (ON(i).gt.thres) then
    !             ONaux(occON) = ON(i)
    !             occON = occON + 1
    !         end if
    !     end do
    !     ! Store these ON in ON red
    !     allocate(ONred(2*occON), stat=ierr)
    !     if (ierr .ne. 0) stop 'Eee_PNOF: Error in allocation of ONred'
    !     do i = 1, occON
    !         ONred(i) = ONaux(i)/2.d0
    !         ONred(2*occON-i+1) = ONaux(i)/2.d0
    !     end do
    !
    !     write(*,*) 'Number of occupied orbitals:', occON
    !     write(*,*) 'ONred:'
    !     write(*,*) ONred
    !     write(*,*) 'sum(ON):', sum(ONred)
    !
    !     ! Divide D' into N/2 pairs fulfilling
    !     ! np + nq = 1 with p,q coupled orbitals in a pair P
    !     ! for an element of D' find another element such that its sum = 1
    !     ! and store them in a pair P
    !     ! remove those orbitals and loop over D' until finding all P pairs
    !     allocate(pairs(2,occON), stat=ierr)
    !     if (ierr .ne. 0) stop 'Eee_PNOF: Error in allocation of pairs'
    !     pairs = 0.d0
    !     ! each pair P is stored as a column in 'pairs'
    !     P = 1
    !     pairs(1,P) = ONred(P)
    !     pair_indeces = 
    !     do i = 2, 2*occON
    !         if (  (ONred(i) + pairs(i,P) > 1.d0 - thres) .or. &
    !             & (ONred(i) + pairs(i,P) < 1.d0 + thres)) then
    !             pairs(2,P) = ONred(i)
    !             P = P + 1
    !         end if
    !     end do
    !     
    !
    !     ! Compute EeePNOF5
    !     !
    !     return
    ! end subroutine Eee_PNOF 

    subroutine energy(E1, E2, Etot)
        implicit none
        real(kind=8), intent(in) :: E1, E2
        real(kind=8), intent(out) :: Etot
        !
        Etot = E1 + E2
        !
        return
    end subroutine energy 

end module density_matrices 
