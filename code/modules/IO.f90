! +--------------------------------------------+
! | Author: Jose Antonio Quinonero Gris        |
! | Creation date: Tuesday 17:27:38 30/04/2024 |
! +--------------------------------------------+
module io 
    ! 
    ! Input/Output module
    !
    use declarations
    use strings
    use matmod
    !
    implicit none
    !
    contains

    subroutine read_input(ifn, d1fn, d2fn, h1fn, h2fn)
        !
        ! ifn: input file name
        ! d1fn: one-electron density matrix file name
        ! d2fn: two-electron density matrix file name
        ! h1fn: one-electron integrals file name
        ! h2fn: two-electron integrals file name
        !
        implicit none
        character(len=*), intent(in) :: ifn
        character(len=80), intent(out) :: d1fn, d2fn, h1fn, h2fn
        !
        integer :: ios, nlu 
        character(len=80) :: datafolder 

        ! Open input file

        open(newunit=nlu, file=trim(ifn), iostat=ios)
        if (ios /= 0) stop "read_input: Error opening file trim(ifn)"

        ! Read input file
        read(nlu,*)       ! skip line
        read(nlu,*) datafolder
        read(nlu,*)
        read(nlu,*) d1fn
        read(nlu,*)
        read(nlu,*) d2fn
        read(nlu,*)
        read(nlu,*) h1fn
        read(nlu,*)
        read(nlu,*) h2fn

        ! Add full path to the filenames
        d1fn = trim(datafolder) // '/' // trim(d1fn)
        d2fn = trim(datafolder) // '/' // trim(d2fn)
        h1fn = trim(datafolder) // '/' // trim(h1fn)
        h2fn = trim(datafolder) // '/' // trim(h2fn)

        ! Close file
        close(unit=nlu, iostat=ios, status="keep")
        if (ios /= 0) stop "read_input: Error closing file unit nlu"
        !
        return
    end subroutine read_input 

    subroutine read_2D_matrix(fn, is_triangular, Amat)
        !
        ! fn: matrix input file name
        ! Amat: any 2D matrix
        !
        implicit none
        character(len=*), intent(in) :: fn
        logical :: is_triangular
        real(kind=8), dimension(:,:), allocatable, intent(out) :: Amat
        !
        integer :: i, j, ios, nlu, n, iind, jind
        real(kind=8) :: xint 

        ! Open file
        open(newunit=nlu, file=trim(fn), iostat=ios)
        if (ios /= 0) stop "read_2D_matrix: Error opening file trim(fn)"

        ! Read dimensions of Amat
        read(nlu,*) 
        read(nlu,*) n
        read(nlu,*) 

        ! Allocate Amat
        allocate(Amat(n,n), stat=ios)
        if (ios /= 0) stop "read_2D_matrix: Error in allocation of Amat"

        ! Read the triangular part of Amat
        iind = -1
        jind = -1
        do
            read(nlu,*,iostat=ios) jind, iind, xint
            if (ios /= 0) stop 'read_2D_matrix: Error reading file.'
            if (iind.eq.0 .or. jind.eq.0) exit
            Amat(iind,jind) = xint
        end do

        ! Mirror the upper triangular
        if (is_triangular .eqv. .true.) then
            do j=1, n-1
                do i=j+1, n
                    Amat(i,j) = Amat(j,i)
                end do
            end do
        end if

        ! Close file
        close(unit=nlu, iostat=ios, status="keep")
        if (ios /= 0) stop "read_2D_matrix: Error closing file unit nlu"
        
        !
        return
    end subroutine read_2D_matrix 

    subroutine read_4D_matrix(fn, is_triangular, Amat)
        !
        ! fn: matrix input file name
        ! Amat: any 4D matrix
        !
        implicit none
        character(len=*), intent(in) :: fn
        logical :: is_triangular
        real(kind=8), dimension(:,:,:,:), allocatable, intent(out) :: Amat
        !
        integer :: i, j, k, l, ios, nlu, n, iind, jind, kind, lind
        real(kind=8) :: xint 

        ! Open file
        open(newunit=nlu, file=trim(fn), iostat=ios)
        if (ios /= 0) stop "read_4D_matrix: Error opening file trim(fn)"

        ! Read dimensions of Amat
        read(nlu,*) 
        read(nlu,*) n
        read(nlu,*) 

        ! Allocate Amat
        allocate(Amat(n,n,n,n), stat=ios)
        if (ios /= 0) stop "read_4D_matrix: Error in allocation of Amat"
        Amat = 0.d0

        ! Read the triangular part of Amat
        iind = -1
        jind = -1
        kind = -1
        lind = -1
        do
            read(nlu,*,iostat=ios) kind, lind, iind, jind, xint
            if (ios /= 0) stop 'read_4D_matrix: Error reading file.'
            if (kind.eq.0 .or. lind.eq.0 .or. iind.eq.0 .or. jind.eq.0) exit

            ! Note: the integrals are written from Dalton in chemistry notation
            ! xint = (ab|cd) -> <ac|bd>
            ! in this program, they are stored in braket notation
            ! TODO
            Amat(iind,jind,kind,lind) = xint
            ! Amat(iind,kind,jind,lind) = xint
        end do

        ! Simetrize the 4D tensor
        if (is_triangular .eqv. .true.) then
            do k=1, n
                do l=k, n
                    do i=1, n
                        do j=i, n
                            Amat(j,i,l,k) = Amat(i,j,k,l)
                            Amat(i,j,l,k) = Amat(i,j,k,l)
                            Amat(j,i,k,l) = Amat(i,j,k,l)
                        end do
                    end do
                end do
            end do
        end if

        ! Close file
        close(unit=nlu, iostat=ios, status="keep")
        if (ios /= 0) stop "read_4D_matrix: Error closing file unit nlu"
        
        !
        return
    end subroutine read_4D_matrix 

    subroutine read_neptunus(fn, is_triangular, H1mat, H2mat)
        !
        ! fn: matrix input file name
        ! H2mat: any 4D matrix
        !
        implicit none
        character(len=*), intent(in) :: fn
        logical :: is_triangular
        real(kind=8), dimension(:,:), allocatable, intent(out) :: H1mat
        real(kind=8), dimension(:,:,:,:), allocatable, intent(out) :: H2mat
        !
        logical :: found_norb
        integer :: i, j, k, l, pos, ios, nlu, n, iind, jind, kind, lind
        real(kind=8) :: xint 
        character(len=100) :: line 

        ! Open file
        open(newunit=nlu, file=trim(fn), iostat=ios)
        if (ios /= 0) stop "read_4D_matrix: Error opening file trim(fn)"

        ! Read dimensions of H2mat
        found_norb = .false.

        ! Read the file
        read(nlu, *)  ! skip first line
        read(nlu, '(A)', iostat=ios) line
        if (ios /= 0) stop 'read_4D_matrix_neptunus: Error reading file.'

        ! Search for "NORB="
        pos = index(line, 'NORB=')
        if (pos /= 0) then
            read(line(pos+5:), *) n  ! Read NORB
        end if

        ! Skip following 3 lines
        do i=1,3; read(nlu,*) ; end do

        ! Allocate H2mat
        allocate(H1mat(n,n), stat=ios)
        if (ios /= 0) stop "read_neptunus: Error in allocation of H1mat"
        H1mat = 0.d0

        ! Allocate H2mat
        allocate(H2mat(n,n,n,n), stat=ios)
        if (ios /= 0) stop "read_neptunus: Error in allocation of H2mat"
        H2mat = 0.d0

        ! Read the triangular part of H2mat
        iind = -1
        jind = -1
        kind = -1
        lind = -1
        do
            read(nlu,*,iostat=ios) xint, kind, lind, iind, jind
            if (ios /= 0) stop 'read_neptunus: Error reading file.'

            if (kind.eq.0 .and. lind.eq.0 .and. iind.eq.0 .and. jind.eq.0) exit  ! end of file

            if (iind.ne.0 .and. jind.ne.0) then
                H2mat(iind,jind,kind,lind) = xint
            else
                ! H1 mat when the last 2 indeces are 0
                H1mat(kind,lind) = xint
            end if
        end do

        if (is_triangular .eqv. .true.) then
            do k=1, n
                do l=k, n
                    ! Simetrize the triangular H1mat matrix
                    H1mat(l,k) = H1mat(k,l)
                    do i=1, n
                        do j=i, n
                            ! Simetrize the 4D tensor H2mat
                            H2mat(j,i,l,k) = H2mat(i,j,k,l)
                            H2mat(i,j,l,k) = H2mat(i,j,k,l)
                            H2mat(j,i,k,l) = H2mat(i,j,k,l)
                        end do
                    end do
                end do
            end do
        end if

        ! Close file
        close(unit=nlu, iostat=ios, status="keep")
        if (ios /= 0) stop "read_neptunus: Error closing file unit nlu"
        
        !
        return
    end subroutine read_neptunus

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

    subroutine print_energies
        integer :: eouf 
        !
        open(newunit=eouf, file="energies.dat")

        write(eouf,'(a)') 'ED1D2      ELS        BBC1       BBC2       BBC3       BBC3M'
        write(eouf,'(*(f9.4,2x))') E_exact, EELS, EBBC1, EBBC2, EBBC3, EBBC3M

        close(eouf)
        !
        return
    end subroutine print_energies 

end module io 
