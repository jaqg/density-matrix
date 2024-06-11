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
    !
    implicit none
    !
    contains
    
    function trace(Amat) result(res)
        !
        ! Function to compute the trace of a matrix, Amat
        !
        implicit none
        real(kind=8), dimension(:,:), intent(in) :: Amat
        real(kind=8) :: res 
        !
        integer :: i 
    
        res = 0.d0 
        do i = 1, size(Amat, dim=1)
            res = res + Amat(i,i)
        end do 
    
    end function trace 

    subroutine diag_2D_mat(Amat, EVvec)
        !
        ! Subroutine to diagonalize a 2D matrix using the LAPACK routine
        ! DSYEV()
        !
        implicit none
        real(kind=8), dimension(:,:), intent(in) :: Amat
        real(kind=8), dimension(:), allocatable, intent(out) :: EVvec
        !
        real(kind=8), dimension(:,:), allocatable :: Aaux
        real(kind=8), dimension(:), allocatable :: work
        integer :: n, ierr, lwork, info

        ! Save the original matrix
        allocate(Aaux(size(Amat, dim=1),size(Amat, dim=2)), stat=ierr)
        if (ierr .ne. 0) stop 'diag_2D_mat: Error in allocation of Aaux'
        !
        Aaux = Amat        

        ! Check for inconsistency in the dimensions
        n = size(Aaux, dim=1)
        if (n .ne. size(Aaux, dim=2)) &
        & stop 'diag_2D_mat: Error: Inconsistent dimensions'

        ! Set the size of the WORK vector for the diagonalization algorithm
        ! LWORK=3*N-1  ! minimum
        lwork = 3 * n

        ! Allocate arrays
        allocate(work(lwork), stat=ierr)
        if (ierr .ne. 0) stop 'diag_2D_mat: Error in allocation of work'
        !
        allocate(EVvec(n), stat=ierr)
        if (ierr .ne. 0) stop 'diag_2D_mat: Error in allocation of EVvec'
        
        ! Call the diagonalization subroutine
        ! checkear dsyevd
        CALL DSYEV( 'V', 'U', n, Aaux, n, EVvec, work, lwork, info )
        if (info.lt.0) then 
            write(*,'(a,i0,a)') 'diag_2D_mat: Error: the ',info, &
            & '-th argument had an illegal value'
            stop
        else if (info.gt.0) then
            write(*,'(a,i0,a)') 'diag_2D_mat: Error: the algorithm failed to &
            & converge; ', info, ' off-diagonal elements of an intermediate &
            & tridiagonal form did not converge to zero.'
            stop
        end if
        
        !
        return
    end subroutine diag_2D_mat 

    subroutine MO_to_SO(AMO, ASO)
        implicit none
        real(kind=8), dimension(:,:), intent(in) :: AMO
        real(kind=8), dimension(:,:), allocatable, intent(out) :: ASO
        !
        integer :: ierr, i, j, n

        ! Get size of the matrix & asure it is square
        n = size(AMO, dim=1)
        if (n .ne. size(AMO, dim=2)) &
        & stop 'MO_to_SO: Error: input matrix not square'
        
        ! Allocate ASO
        allocate(ASO(2*n, 2*n), stat=ierr)
        if (ierr .ne. 0) stop 'MO_to_SO: Error in allocation of ASO'
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
    end subroutine MO_to_SO 

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
            !
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

    subroutine Eee_PNOF(ON, thres)
        implicit none
        real(kind=8), dimension(:), intent(in) :: ON
        real(kind=8), intent(in) :: thres 
        !
        integer :: ierr, i, j, n, occON
        real(kind=8), dimension(:), allocatable :: ONaux, ONred
        real(kind=8), dimension(:,:), allocatable :: pairs

        !
        ! For PNOF5:
        !
        ! Extract the density matrix only with occupied orbitals, i.e.
        ! D'_ij <- D^(1)_ij > 0
        ! which, for the ON vector read
        ! ON'_i <- ON_i > 0
        ! Get dimensions of D1
        n = size(ON)
        allocate(ONaux(n), stat=ierr)
        if (ierr .ne. 0) stop 'Eee_PNOF: Error in allocation of ONaux'
        ONaux = 0.d0
        occON = 1
        do i = 1, n
            if (ON(i).gt.thres) then
                ONaux(occON) = ON(i)
                occON = occON + 1
            end if
        end do
        ! Store these ON in ON red
        allocate(ONred(2*occON), stat=ierr)
        if (ierr .ne. 0) stop 'Eee_PNOF: Error in allocation of ONred'
        do i = 1, occON
            ONred(i) = ONaux(i)/2.d0
            ONred(2*occON-i+1) = ONaux(i)/2.d0
        end do

        write(*,*) 'Number of occupied orbitals:', occON
        write(*,*) 'ONred:'
        write(*,*) ONred
        write(*,*) 'sum(ON):', sum(ONred)

        ! Divide D' into N/2 pairs fulfilling
        ! np + nq = 1 with p,q coupled orbitals in a pair P
        ! for an element of D' find another element such that its sum = 1
        ! and store them in a pair P
        ! remove those orbitals and loop over D' until finding all P pairs
        allocate(pairs(2,occON), stat=ierr)
        if (ierr .ne. 0) stop 'Eee_PNOF: Error in allocation of pairs'
        pairs = 0.d0
        ! each pair P is stored as a column in 'pairs'
        P = 1
        pairs(1,P) = ONred(P)
        pair_indeces = 
        do i = 2, 2*occON
            if (  (ONred(i) + pairs(i,P) > 1.d0 - thres) .or. &
                & (ONred(i) + pairs(i,P) < 1.d0 + thres)) then
                pairs(2,P) = ONred(i)
                P = P + 1
            end if
        end do
        

        ! Compute EeePNOF5
        !
        return
    end subroutine Eee_PNOF 

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
