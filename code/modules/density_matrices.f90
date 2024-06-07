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
        do i = 1, n
            do j = 1, n
                E1 = E1 + H1(i,j) * D1(i,j)
            end do
        end do

        !
        return
    end subroutine H1D1_energy 

    ! subroutine H1D1_energy(D1, H1, E1)
    !     !
    !     ! Compute the non-ee part of the energy: tr[D^(1) H^(1)]
    !     !
    !     ! D1: one-electron density matrix
    !     ! H1: one-electron integral matrix
    !     !
    !     ! E1: energy
    !     !
    !     implicit none
    !     real(kind=8), dimension(:,:), intent(in) :: D1, H1
    !     real(kind=8), intent(out) :: E1
    !     !
    !     integer :: i, j, n 
    !     real(kind=8) :: sum1, sum2 
    !
    !     ! Get the dimensions
    !     n = size(D1,1)
    !     if (n .ne. size(H1,1) .or. n .ne. size(D2,1) .or. n .ne. size(H2,1)) &
    !     & stop 'HF_energy: Error: Inconsistent dimensions'
    !
    !     ! Compute the first sum
    !     sum1 =  0.d0
    !     do i = 1, n
    !         do j = 1, n
    !             sum1 = sum1 + H1(i,j) * D1(i,j)
    !         end do
    !     end do
    !
    !     ! Compute the second sum
    !     sum2 =  0.d0
    !     do i = 1, n
    !         do j = 1, n
    !             do k = 1, n
    !                 do l = 1, n
    !                     sum2 = sum2 + H2(i,j,k,l) * D2(i,j,k,l)
    !                 end do
    !             end do
    !         end do
    !     end do
    !
    !     ! Compute the energy
    !     EHF = sum1 + 0.5d0 * sum2
    !     
    !     !
    !     return
    ! end subroutine HF_energy 
end module density_matrices 
