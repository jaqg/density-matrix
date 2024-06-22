module PNOF
    contains
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

end module PNOF
