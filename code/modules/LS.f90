module LS
    use density_matrices

    implicit none
    contains

    function GLS(p,q,np,nq) result(res)
        implicit none
        integer(kind=8), intent(in) :: p,q
        real(kind=8), intent(in) :: np,nq
        real(kind=8) :: res 

        if (p.eq.q) then 
            res = np
        elseif ((p.eq.1 .and. q.gt.1) .or. (p.gt.1 .and. q.eq.1)) then 
            res = - dsqrt(np * nq)
        else
            res = dsqrt(np * nq)
        end if

    end function GLS 

    subroutine calc_Eee_LS(ONin, H2, nisht, nasht, eina, eact, ecross, Eee)
        implicit none
        real(kind=8), dimension(:), intent(in) :: ONin
        real(kind=8), dimension(:,:,:,:), intent(in) :: H2
        integer(kind=8), intent(in) :: nisht, nasht 
        real(kind=8), intent(out) :: eina, eact, ecross, Eee
        !
        integer(kind=8) :: i, j, u, uu, v, vv
        integer :: n , ierr 
        real(kind=8) :: ni, nj, nu, nv
        real(kind=8), dimension(:), allocatable :: ON

        ! Check for inconsistencies in dimension
        n = size(ONin)
        if (n.ne.size(H2, dim=1) .or. &
          & n.ne.size(H2, dim=2) .or. &
          & n.ne.size(H2, dim=3) .or. &
          & n.ne.size(H2, dim=4)) &
          stop 'calc_Eee_LS: inconsistent dimensions'

        ! Check that ON is in spin-orbital basis
        ! call check_SO_ON(ONin)

        ! Check that ON is sorted in descending order
        allocate(ON(n), stat=ierr)
        if (ierr .ne. 0) stop 'calc_Eee_LS: Error in allocation of ON'
        ON = ONin
        call asure_descending(ONin, ON)
        ! TODO: implementar en caso de que sea false para que no se pierda
        ! la correspondencia en el indexado entre ocupaciones-integrales

        !
        ! Inactive energy
        !
        eina = 0.d0
        do i = 1, nisht
            do j = 1, nisht
                ni = ON(i) ; nj = ON(j)
                eina = eina + 0.5d0 * GLS(i,j,ni,nj) * H2(i,j,j,i)
            end do
        end do

        !
        ! Inactive-active energy
        !
        ecross = 0.d0
        do i = 1, nisht
            do u = 1, nasht
                uu = nisht + u
                ni = ON(i) ; nu = ON(uu)
                ecross = ecross + 0.5d0 * GLS(i,uu,ni,nu) * H2(i,uu,uu,i)
            end do
        end do
        do u = 1, nasht
            uu = nisht + u
            do i = 1, nisht
                ni = ON(i) ; nu = ON(uu)
                ecross = ecross + 0.5d0 * GLS(uu,i,nu,ni) * H2(uu,i,i,uu)
            end do
        end do
        !
        ! Active-active energy
        !
        eact = 0.d0
        do u = 1, nasht
            uu = nisht + u
            do v = 1, nasht
                vv = nisht + v
                nu = ON(uu) ; nv = ON(vv)
                eact = eact + 0.5d0 * GLS(uu,vv,nu,nv) * H2(uu,vv,vv,uu)
            end do
        end do

        ! Compute Eee
        Eee = eina + ecross + eact

    end subroutine calc_Eee_LS 

    subroutine Eee_LS(ONin, H2, Eee)
        implicit none
        real(kind=8), dimension(:), intent(in) :: ONin
        real(kind=8), dimension(:,:,:,:), intent(in) :: H2
        real(kind=8), intent(out) :: Eee
        !
        integer(kind=8) :: i, j, n
        integer :: ierr 
        real(kind=8) :: sum1, ni, nj
        real(kind=8), dimension(:), allocatable :: ON

        ! Check for inconsistencies in dimension
        n = size(ONin)
        if (n.ne.size(H2, dim=1) .or. &
          & n.ne.size(H2, dim=2) .or. &
          & n.ne.size(H2, dim=3) .or. &
          & n.ne.size(H2, dim=4)) &
          stop 'Eee_LS: Error: inconsistent dimensions'

        ! Check that ON is in spin-orbital basis
        call check_SO_ON(ON)

        ! Check that ON is sorted in descending order
        allocate(ON(n), stat=ierr)
        if (ierr .ne. 0) stop 'Eee_LS: Error in allocation of ON'
        ON = ONin
        call asure_descending(ONin, ON)

        ! Compute the sum
        sum1 = 0.d0
        do i = 1, n
            do j = 1, n
                ni = ON(i) ; nj = ON(j)
                sum1 = sum1 + GLS(i,j,ni,nj) * H2(i,j,j,i)
            end do
        end do

        ! Compute Eee
        Eee = 0.5d0 * sum1

    end subroutine Eee_LS 

    subroutine D2_LS(ON, D2)
        !
        ! Subroutine to compute the LS 2-RDM 
        !
        implicit none
        real(kind=8), dimension(:), intent(in) :: ON
        real(kind=8), dimension(:,:,:,:), allocatable, intent(out) :: D2
        !
        integer(kind=8) :: p, q, n
        integer :: ierr
        real(kind=8) :: np, nq 

        n = size(ON)
        allocate(D2(n,n,n,n), stat=ierr)
        if (ierr .ne. 0) stop 'D2_LS: Error in allocation of D2'
        D2 = 0.d0

        do q = 1, n
            do p = 1, n
                np = ON(p) ; nq = ON(q)
                D2(p,q,q,p) = GLS(p,q,np,nq)
            end do
        end do

        D2 = 0.5d0 * D2
        
        !
        return
    end subroutine D2_LS 

end module LS
