module BBC
    use density_matrices

    implicit none
    contains

    real(kind=8) function kronecker(i,j)
        implicit none
        integer, intent(in) :: i,j
    
        if (i .eq. j) then
            kronecker = 1.d0
        else
            kronecker = 0.d0
        end if
    
    end function kronecker 


    function is_weak(np) result(isweak)
        !
        ! Function to check if p is a weakly occupied orbital based on np
        !
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

    function is_bonding(np) result(isbonding)
        !
        ! Function to check if p is a bonding orbital based on np
        !
        implicit none
        real(kind=8), intent(in) :: np
        logical :: isbonding
        !
        isbonding = .false.
        if (np.gt.0.9d0 .and. np.lt.1.d0) then
            isbonding = .true.
        end if
    endfunction is_bonding

    function is_antibond(np) result(isantibond)
        !
        ! Function to check if p is an antibonding orbital based on np
        !
        implicit none
        real(kind=8), intent(in) :: np
        logical :: isantibond
        !
        isantibond = .false.
        if (np.gt.0.d0 .and. np.lt.0.1d0) then
            isantibond = .true.
        end if
    endfunction is_antibond

    function is_frontier(np) result(isfrontier)
        !
        ! Function to check if p is a frontier orbital based on np
        !
        implicit none
        real(kind=8), intent(in) :: np
        logical :: isfrontier

        isfrontier = .false.
        if (is_bonding(np) .or. is_antibond(np)) then
            isfrontier = .true.
        end if
    endfunction is_frontier

    function is_degenerate(np,nq) result(isdeg)
        !
        ! Function to check if p,q are degenerate orbitals based on np,nq
        !
        implicit none
        real(kind=8), intent(in) :: np, nq
        logical :: isdeg

        isdeg = .false.
        if (dabs(np-nq).lt.1.d-5) then
            isdeg = .true.
        end if
    endfunction is_degenerate

    function in_subset(p, subset_ind) result(insub)
        !
        ! Function to check if p index in is subset_ind
        !
        implicit none
        integer, intent(in) :: p
        integer, dimension(:), intent(in) :: subset_ind
        logical :: insub
        !
        integer :: i, n 

        n = size(subset_ind)

        insub = .false.
        do i = 1, n
            if (p .eq. subset_ind(i)) then
                insub = .true.
            end if
        end do
    endfunction in_subset

    subroutine add_to_set(set, set_ind, element, index)
        !
        ! Subroutine to add the element 'element' to the end of 'set', and
        ! the index 'index' to the end of 'set_ind'.
        !
        implicit none
        real(kind=8), allocatable, intent(inout) :: set(:)
        integer, allocatable, intent(inout) :: set_ind(:)
        real(kind=8), intent(in) :: element
        integer, intent(in) :: index
        !
        integer :: n
        real(kind=8), dimension(:), allocatable :: set_old
        integer, dimension(:), allocatable :: set_ind_old

        n = size(set)
        allocate(set_old(n), set_ind_old(n))
        set_old = set
        set_ind_old = set_ind
        if (n == 0) then
            deallocate(set, set_ind)
            allocate(set(1))
            allocate(set_ind(1))
            set(1) = element
            set_ind(1) = index
        else
            deallocate(set, set_ind)
            allocate(set(n+1))
            allocate(set_ind(n+1))
            set(1:n) = set_old
            set(n+1) = element
            set_ind(1:n) = set_ind_old
            set_ind(n+1) = index
        end if
    end subroutine add_to_set

    subroutine subdivide_S(On, Sb, Sb_ind, Sc, Sc_ind)
        !
        ! Subroutine to apply the modification by Lathiotakis and Marques:
        ! Subdivide S into:
        ! - Sb for the degenerate bonding orbitals 
        ! - Sc for the rest
        !
        implicit none
        real(kind=8), dimension(:), intent(in) :: On
        real(kind=8), dimension(:), allocatable, intent(out) :: Sb, Sc
        integer, dimension(:), allocatable, intent(out) :: Sb_ind, Sc_ind
        !
        integer :: i, j, n, ierr
        logical :: is_deg
        real(kind=8) :: ni, nj 

        ! Get the size of the input set
        n = size(On)

        ! Initialize subsets and index arrays as empty
        allocate(Sb(0), stat=ierr)
        if (ierr .ne. 0) stop 'subdivide_S: Error in allocation of Sb'
        !
        allocate(Sc(0), stat=ierr)
        if (ierr .ne. 0) stop 'subdivide_S: Error in allocation of Sc'
        !
        allocate(Sb_ind(0), stat=ierr)
        if (ierr .ne. 0) stop 'subdivide_S: Error in allocation of Sb_ind'
        !
        allocate(Sc_ind(0), stat=ierr)
        if (ierr .ne. 0) stop 'subdivide_S: Error in allocation of Sc_ind'

        ! Iterate over the orbitals in ON
        do i = 1, n
            if (in_subset(i,Sb_ind)) cycle
            ni = ON(i)
            if (.not.is_weak(ni) .and. is_bonding(ni) .and. i.lt.n) then
                ! Check for degeneracy
                is_deg = .false.
                do j = i+1, n
                    nj = ON(j)
                    if (.not.is_weak(nj) .and. &
                        &is_bonding(nj) .and. is_degenerate(ni, nj)) then
                        is_deg = .true.
                        exit
                    end if
                end do
                if (is_deg) then
                    ! Add ni to Sb if it's bonding and degenerate with nj
                    call add_to_set(Sb, Sb_ind, ni, i)
                    ! Add nj to Sb if it's bonding and degenerate with ni
                    call add_to_set(Sb, Sb_ind, nj, j)
                else
                    ! Add to Sc otherwise
                    call add_to_set(Sc, Sc_ind, ni, i)
                end if
            else
                ! Add to Sc if it's not a bonding orbital
                call add_to_set(Sc, Sc_ind, ni, i)
            end if
        end do
        return
    end subroutine subdivide_S 

    subroutine subdivide_W(ON, Wa, Wa_ind, Wh, Wh_ind)
        !
        ! Subroutine to apply the modification by Lathiotakis and Marques:
        ! Subdivide W into:
        ! - Wa for the degenerate antibonding orbitals 
        ! - Wh for the rest
        !
        implicit none
        real(kind=8), dimension(:), intent(in) :: On
        real(kind=8), dimension(:), allocatable, intent(out) :: Wa, Wh
        integer, dimension(:), allocatable, intent(out) :: Wa_ind, Wh_ind
        !
        integer :: i, j, n, ierr
        logical :: is_deg
        real(kind=8) :: ni, nj 

        ! Get the size of the input set
        n = size(On)

        ! Initialize subsets and index arrays as empty
        allocate(Wa(0), stat=ierr)
        if (ierr .ne. 0) stop 'subdivide_S: Error in allocation of Wa'
        !
        allocate(Wh(0), stat=ierr)
        if (ierr .ne. 0) stop 'subdivide_S: Error in allocation of Wh'
        !
        allocate(Wa_ind(0), stat=ierr)
        if (ierr .ne. 0) stop 'subdivide_S: Error in allocation of Wa_ind'
        !
        allocate(Wh_ind(0), stat=ierr)
        if (ierr .ne. 0) stop 'subdivide_S: Error in allocation of Wh_ind'

        ! Iterate over the orbitals in ON
        do i = 1, n
            if (in_subset(i,Wa_ind)) cycle
            ni = ON(i)
            if (is_weak(ni) .and. is_antibond(ni) .and. i.lt.n) then
                ! Check for degeneracy
                is_deg = .false.
                do j = i+1, n
                    nj = ON(j)
                    if (is_weak(nj) .and. &
                        &is_antibond(nj) .and. is_degenerate(ni, nj)) then
                        is_deg = .true.
                        exit
                    end if
                end do
                if (is_deg) then
                    ! Add ni to Wa if it's bonding and degenerate with nj
                    call add_to_set(Wa, Wa_ind, ni, i)
                    ! Add nj to Wa if it's bonding and degenerate with ni
                    call add_to_set(Wa, Wa_ind, nj, j)
                else
                    ! Add to Wh otherwise
                    call add_to_set(Wh, Wh_ind, ni, i)
                end if
            else
                ! Add to Wh if it's not a bonding orbital
                call add_to_set(Wh, Wh_ind, ni, i)
            end if
        end do
        return
    end subroutine subdivide_W

    subroutine calc_GBBC(level, ON, GBBC)
        implicit none
        character(len=*), intent(in) :: level 
        real(kind=8), dimension(:), intent(in) :: ON
        real(kind=8), dimension(:,:), allocatable, intent(out) :: GBBC 
        !
        integer :: i, j, n, ierr
        real(kind=8) :: ni, nj 
        real(kind=8), dimension(:), allocatable :: Sb, Sc, Wa, Wh 
        integer, dimension(:), allocatable :: Sb_ind, Sc_ind, Wa_ind, Wh_ind 

        n = size(ON)

        ! Compute the GBBC matrix
        allocate(GBBC(n,n), stat=ierr)
        if (ierr .ne. 0) stop 'calc_BBC: Error in allocation of GBBC'
        GBBC = 0.d0

        if (level.eq.'1' .or. level.eq.'BBC1') then
            ! Corrected BB functional, BBC1
            do i = 1, n
                do j = 1, n
                    ni = ON(i) ; nj = ON(j)
                    if (j.ne.i .and. is_weak(ni) .and. is_weak(nj)) then
                        ! For p/=q and p,q are weakly occupied (np,nq < 1/2)
                        GBBC(i,j) = dsqrt(ni * nj)
                    else
                        GBBC(i,j) = - dsqrt(ni * nj)
                    end if
                end do 
            end do
        elseif (level.eq.'2' .or. level.eq.'BBC2') then
            do i = 1, n
                do j = 1, n
                    ni = ON(i) ; nj = ON(j)
                    if (j.ne.i .and. is_weak(ni) .and. is_weak(nj)) then
                        ! For p/=q and p,q are weakly occupied (np,nq < 1/2)
                        GBBC(i,j) = dsqrt(ni * nj)
                    elseif (j.ne.i .and. .not.is_weak(ni) &
                        &          .and. .not.is_weak(nj)) then
                        ! For p/=q and p,q are strongly occupied (np,nq > 1/2)
                        GBBC(i,j) = -ni * nj
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
                        &   (is_weak(ni) .and. is_frontier(nj) &
                        &                .and.is_weak(nj)) .or. &
                        ! Or p weak, q frontier (weak) 
                        &   (is_frontier(ni) .and. is_weak(ni) &
                        &                    .and. is_weak(nj))&
                        &) then
                        ! Or p frontier (weak), q weak 
                        GBBC(i, j) = dsqrt(ni * nj)
                    elseif ((i.ne.j .and. .not.is_weak(ni) .and. &
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
        elseif (level.eq.'3M' .or. level.eq.'BBC3M') then
            !
            ! Modified BBC3 by Lathiotakis and Marques
            !
            ! Subdivide S into Sb,Sc and W into Wa, Wh
            call subdivide_S(ON, Sb, Sb_ind, Sc, Sc_ind)
            call subdivide_W(ON, Wa, Wa_ind, Wh, Wh_ind)

            do i = 1, n
                do j = 1, n
                    ni = ON(i) ; nj = ON(j)
                    if (j.ne.i .and. is_weak(ni) .and. is_weak(nj)) then
                        ! For p/=q and p,q are weakly occupied (np,nq < 1/2)
                        GBBC(i,j) = dsqrt(ni * nj)
                    elseif (j.ne.i .and. .not.is_weak(ni) &
                        &          .and. .not.is_weak(nj)) then
                        ! For p/=q and p,q are strongly occupied (np,nq > 1/2)
                        GBBC(i,j) = -ni * nj
                    elseif (&
                        &(in_subset(i,Sc_ind) .and. in_subset(j,Wa_ind))&
                        &.or.&
                        &(in_subset(i,Wa_ind) .and. in_subset(j,Sc_ind))&
                        &) then
                        ! For p in Sc and q in Wa, or p in Wa and q in Sc
                        GBBC(i,j) = -ni * nj
                    elseif (j.ne.i .and. &
                        &(in_subset(i,Sc_ind) .or. in_subset(i,Wh_ind))&
                        &.and.&
                        &(in_subset(j,Sc_ind) .or. in_subset(j,Wh_ind))&
                        &) then
                        ! For p/= q and p,q in {Sc + Wh}
                        GBBC(i,j) = - ni**2
                    else
                        GBBC(i,j) = - dsqrt(ni * nj)
                    end if
                end do 
            end do
        else
            stop 'calc_BBC: Error: Invalid level. Available: BBC1, BBC2, BBC3, BBC3M'
        end if
    end subroutine calc_GBBC 
    
    subroutine Eee_BBC(level, ON, H2, Eee)
        !
        ! Subroutine to compute the Eee for the BBCn approximations in
        ! the spin-orbital basis with
        ! 0 <= np <= 1 forall p
        !
        implicit none
        character(len=*), intent(in) :: level 
        real(kind=8), dimension(:), intent(in) :: ON
        real(kind=8), dimension(:,:,:,:), intent(in) :: H2
        real(kind=8), intent(out) :: Eee
        !
        integer :: i, j, n
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
        call check_SO_ON(ON)

        ! Compute the GBBC matrix
        call calc_GBBC(level, ON, GBBC)

        ! Compute first sum
        sum1 = 0.d0
        do i = 1, n
            do j = 1, n
                ni = ON(i) ; nj = ON(j)
                sum1 = sum1 + ni * nj * H2(i,i,j,j)
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

        !
    end subroutine Eee_BBC 

    subroutine calc_Eee_BBC(level,ON,H2,nisht,nasht,eina,eact,ecross,Eee)
        !
        ! Subroutine to compute the Eee for the BBCn approximations in
        ! the spin-orbital basis with
        ! 0 <= np <= 1 forall p
        !
        implicit none
        character(len=*), intent(in) :: level 
        real(kind=8), dimension(:), intent(in) :: ON
        real(kind=8), dimension(:,:,:,:), intent(in) :: H2
        integer(kind=8), intent(in) :: nisht, nasht 
        real(kind=8), intent(out) :: eina, eact, ecross, Eee
        !
        integer(kind=8) :: n, i, j, u, uu, v, vv
        real(kind=8) :: ni, nj, nu, nv
        real(kind=8), dimension(:,:), allocatable :: GBBC 

        ! Check for inconsistencies in the dimensions
        n = size(ON)
        if (n.ne.size(H2, dim=1) .or. &
            &   n.ne.size(H2, dim=2) .or. &
            &   n.ne.size(H2, dim=3) .or. &
            &   n.ne.size(H2, dim=4)) &
            & stop 'Eee_BBC: Error: Inconsistent dimensions'

        ! Check that ON is in spin-orbital basis
        call check_SO_ON(ON)

        ! Compute the GBBC matrix
        call calc_GBBC(level, ON, GBBC)

        !
        ! Inactive energy
        !
        eina = 0.d0
        ! First sum
        do i = 1, nisht
            do j = 1, nisht
                ni = ON(i) ; nj = ON(j)
                eina = eina + 0.5d0 * ni * nj * H2(i,i,j,j)
            end do
        end do
        ! Second sum
        do i = 1, nisht
            do j = 1, nisht
                ni = ON(i) ; nj = ON(j)
                eina = eina + 0.5d0 * GBBC(i,j) * H2(i,j,j,i)
            end do
        end do

        !
        ! Inactive-active energy
        !
        ecross = 0.d0
        ! First sum
        do i = 1, nisht
            do u = 1, nasht
                uu = nisht + u
                ni = ON(i) ; nu = ON(uu)
                ecross = ecross + 0.5d0 * ni * nu * H2(i,i,uu,uu)
            end do
        end do
        do u = 1, nasht
            uu = nisht + u
            do i = 1, nisht
                ni = ON(i) ; nu = ON(uu)
                ecross = ecross + 0.5d0 * nu * ni * H2(uu,uu,i,i)
            end do
        end do
        ! Second sum
        do i = 1, nisht
            do u = 1, nasht
                uu = nisht + u
                ni = ON(i) ; nu = ON(uu)
                ecross = ecross + 0.5d0 * GBBC(i,uu) * H2(i,uu,uu,i)
            end do
        end do
        do u = 1, nasht
            uu = nisht + u
            do i = 1, nisht
                ni = ON(i) ; nu = ON(uu)
                ecross = ecross + 0.5d0 * GBBC(uu,i) * H2(uu,i,i,uu)
            end do
        end do
        !
        ! Active-active energy
        !
        eact = 0.d0
        ! First sum
        do u = 1, nasht
            uu = nisht + u
            do v = 1, nasht
                vv = nisht + v
                nu = ON(uu) ; nv = ON(vv)
                eact = eact + 0.5d0 * nu * nv * H2(uu,uu,vv,vv)
            end do
        end do
        ! Second sum
        do u = 1, nasht
            uu = nisht + u
            do v = 1, nasht
                vv = nisht + v
                nu = ON(uu) ; nv = ON(vv)
                eact = eact + 0.5d0 * GBBC(uu,vv) * H2(uu,vv,vv,uu)
            end do
        end do

        ! Compute the energy
        Eee = eina + ecross + eact

        !
    end subroutine calc_Eee_BBC 

    subroutine Eee_BBC_spacial(level, ON, H2, nisht, nasht, &
        & eina, eact, ecross, Eee)
        !
        ! Subroutine to compute the Eee for the BBCn approximations in
        ! the spacial orbital basis with
        ! 0 <= np <= 1 forall p
        !
        implicit none
        character(len=*), intent(in) :: level 
        real(kind=8), dimension(:), intent(in) :: ON
        real(kind=8), dimension(:,:,:,:), intent(in) :: H2
        integer(kind=8), intent(in) :: nisht, nasht 
        real(kind=8), intent(out) :: eina, eact, ecross, Eee
        !
        integer(kind=8) :: i, j, u, uu, v, vv
        integer :: n
        real(kind=8) :: ni, nj, nu, nv
        real(kind=8), dimension(:,:), allocatable :: GBBC 

        ! Check for inconsistencies in the dimensions
        n = size(ON)
        if (n.ne.size(H2, dim=1) .or. &
            &   n.ne.size(H2, dim=2) .or. &
            &   n.ne.size(H2, dim=3) .or. &
            &   n.ne.size(H2, dim=4)) &
            & stop 'Eee_BBC: Error: Inconsistent dimensions'

        ! Compute the GBBC matrix
        call calc_GBBC(level, ON, GBBC)

        !
        ! Inactive energy
        !
        eina = 0.d0
        do i = 1, nisht
            ni = ON(i)
            do j = 1, nisht
                nj = ON(j)
                ! Compute first sum
                ! Coulomb term (pp|qq) -> <pq|pq> 
                ! (since H2 is stored in Dirac notation)
                ! eina = eina + 2.d0 * ni * nj * H2(i,j,i,j)  ! Braket
                eina = eina + 2.d0 * ni * nj * H2(i,i,j,j)  ! Dirac

                ! Compute second sum
                ! Exchange term (pq|qp) -> <pq|qp> 
                ! (since H2 is stored in Dirac notation)
                eina = eina - GBBC(i,j) * H2(i,j,j,i)
            end do
        end do

        !
        ! Inactive-active energy
        !
        ecross = 0.d0
        do i = 1, nisht
            ni = ON(i)
            do u = 1, nasht
                uu = nisht + u
                nu = ON(uu)

                ! Coulomb term (pp|qq) -> <pq|pq> 
                ! ecross = ecross + 2.d0 * ni * nu * H2(i,uu,i,uu)  ! Braket
                ecross = ecross + 2.d0 * ni * nu * H2(i,i,uu,uu)  ! Dirac

                ! Exchange term (pq|qp) -> <pq|qp> 
                ecross = ecross - GBBC(i,uu) * H2(i,uu,uu,i)
            end do
        end do
        do u = 1, nasht
            uu = nisht + u
            nu = ON(uu)
            do i = 1, nisht
                ni = ON(i)

                ! Coulomb term (pp|qq) -> <pq|pq> 
                ! ecross = ecross + 2.d0 * nu * ni * H2(uu,i,uu,i)  ! Braket
                ecross = ecross + 2.d0 * nu * ni * H2(u,uu,i,i)  ! Dirac

                ! Exchange term (pq|qp) -> <pq|qp> 
                ecross = ecross - GBBC(uu,i) * H2(uu,i,i,uu)

            end do
        end do

        !
        ! Active-active energy
        !
        eact = 0.d0
        do u = 1, nasht
            uu = nisht + u
            nu = ON(uu)
            do v = 1, nasht
                vv = nisht + v
                nv = ON(vv)

                ! Coulomb term (pp|qq) -> <pq|pq> 
                ! eact = eact + 2.d0 * nu * nv * H2(uu,vv,uu,vv)  ! Braket
                eact = eact + 2.d0 * nu * nv * H2(uu,uu,vv,vv)  ! Dirac

                ! Exchange term (pq|qp) -> <pq|qp> 
                eact = eact - GBBC(uu,vv) * H2(uu,vv,vv,uu)
            end do
        end do

        ! Compute Eee
        Eee = eina + ecross + eact

        !
    end subroutine Eee_BBC_spacial 

    subroutine D2_BBC(level, ON, D2)
        implicit none
        character(len=*), intent(in) :: level
        real(kind=8), dimension(:), intent(in) :: ON
        real(kind=8), dimension(:,:,:,:), allocatable, intent(out) :: D2
        !
        integer :: p, q, r, s, n, ierr
        real(kind=8) :: np, nq, nr
        real(kind=8), dimension(:,:), allocatable :: GBBC

        ! Compute the GBBC matrix
        call calc_GBBC(level, ON, GBBC)

        ! call print_matrix('GBBC', GBBC(1:6,1:6), 6, '(*(f10.4))')

        n = size(ON)
        allocate(D2(n,n,n,n), stat=ierr)
        if (ierr .ne. 0) stop 'D2_BBC: Error in allocation of D2'
        D2 = 0.d0

        ! do q = 1, n
        !     do p = 1, n
        !         np = ON(p) ; nq = ON(q)
        !         if (p.eq.q) then 
        !             D2(p,p,p,p) = np**2 + GBBC(p,p)
        !         else
        !             D2(p,q,q,p) = GBBC(p,q)
        !         end if
        !     end do
        ! end do
        ! D2 = 0.5d0 * D2  ! ???

        do p = 1, n
            np = ON(p)
            do q = 1, n
                nq = ON(q)
                do r = 1, n
                    nr = ON(r)
                    do s = 1, n
                        ! D2(p,q,r,s) = &
                        ! & 4.d0 * np * nr * kronecker(p,q) * kronecker(r,s) -&
                        ! & 2.d0 * GBBC(p,q) * kronecker(q,r) * kronecker(p,s)
                        !
                        ! In Dirac notation
                        ! D2(p,q,r,s) = &
                        ! & np * nr * kronecker(p,q) * kronecker(r,s) +&
                        ! & GBBC(p,q) * kronecker(q,r) * kronecker(p,s)
                        !
                        ! In Mulliken notation
                        D2(p,q,r,s) = &
                        & np * nq * kronecker(p,r) * kronecker(q,s) +&
                        & GBBC(p,q) * kronecker(q,r) * kronecker(p,s)
                    end do
                end do
            end do
        end do
        D2 = 0.5d0 * D2  ! ???
        
        !
        return
    end subroutine D2_BBC 

    subroutine D2_BBC2(ON, D2)
        implicit none
        real(kind=8), dimension(:), intent(out) :: ON
        real(kind=8), dimension(:,:,:,:), allocatable, intent(out) :: D2
        !
        integer :: n, ierr, p, q, r, s
        real(kind=8) :: np, nq, nr, ns 
        real(kind=8), dimension(:,:), allocatable :: GBBC

        ! Compute the GBBC matrix
        call calc_GBBC('BBC2', ON, GBBC)

        n = size(ON)
        allocate(D2(n,n,n,n), stat=ierr)
        if (ierr .ne. 0) stop 'D2_BBC1: Error in allocation of D2'
        D2 = 0.d0

        do p = 1, n
            np = ON(p)
            do q = 1, n
                nq = ON(q)
                do r = 1, n
                    nr = ON(r)
                    do s = 1, n
                        ns = ON(s)
                        ! D2(p,q,r,s) = np * nq * kronecker(p,r) * kronecker(q,s)&
                        ! & -0.5d0 * (np * nq * kronecker(p,s) * kronecker(q,r) +&
                        ! & np * nq * kronecker(q,r) * kronecker(p,s) )

                        D2(p,q,r,s) = np * nr * kronecker(p,q) * kronecker(r,s)&
                        & + GBBC(p,q) * kronecker(q,r) * kronecker(p,s)
                    end do
                end do
            end do
        end do
        
        !
        return
    end subroutine D2_BBC2

end module BBC
