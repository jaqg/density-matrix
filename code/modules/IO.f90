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

        ! Read the triangular part of Amat
        iind = -1
        jind = -1
        kind = -1
        lind = -1
        do
            read(nlu,*,iostat=ios) kind, lind, iind, jind, xint
            if (ios /= 0) stop 'read_4D_matrix: Error reading file.'
            if (kind.eq.0 .or. lind.eq.0 .or. iind.eq.0 .or. jind.eq.0) exit
            Amat(iind,jind,kind,lind) = xint
        end do

        ! Mirror the upper triangular
        if (is_triangular .eqv. .true.) then
            do l=1, n
                do k=1, n
                    do j=1, n-1
                        do i=j+1, n
                            Amat(i,j,k,l) = Amat(j,i,k,l)
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
        !
        integer :: i, j 

        write(lu,'(a)') label
        do i=1, size(A,1)
            write(lu,trim(fmt)) (A(i,j), j=1, size(A,2))
        end do
        !
        return
    end subroutine print_matrix 

end module io 
