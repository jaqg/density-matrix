! +--------------------------------------------+
! | Author: Jose Antonio Quinonero Gris        |
! | Creation date: Tuesday 09:11:38 16/01/2024 |
! +--------------------------------------------+
module strings

    implicit none
    
    contains

    subroutine separator(ouf, width)
        implicit none
        integer, intent(in) :: ouf, width

        write(ouf,'(*(a))') '+', repeat('-', width-2), '+'

        return
    end subroutine separator 

    subroutine center(uf, mainstr, sep, width)
        implicit none
        integer, intent(in) :: uf
        character(len=*), intent(in) :: mainstr, sep 
        integer, intent(in) :: width
        integer :: numsep 

        if (width .lt. len(mainstr)) stop 'center: width < len(mainstr)'

        numsep = (width - len(mainstr))/2

        if (modulo((width - len(mainstr)),2) .ne. 0) then
            write(uf,'(*(a))') repeat(sep, numsep-1)// &
            & ' '//mainstr//' ' &
            & //repeat(sep, numsep)
        else
            write(uf,'(*(a))') repeat(sep, numsep-1)// &
            & ' '//mainstr//' ' &
            & //repeat(sep, numsep-1)
        end if

        return
    end subroutine center 

    subroutine title(uf, titlestr, width)
        implicit none
        integer, intent(in) :: uf, width
        character(len=*), intent(in) :: titlestr
        integer :: numsep 

        if (width .lt. len(titlestr)) stop 'title: width < len(titlestr)'

        numsep = (width - len(titlestr))/2 - 1

        write(uf,'(*(a))') '+', repeat('=', width-2), '+'
        if (modulo((width - len(titlestr)),2) .ne. 0) then
            write(uf,'(*(a))') '|'//repeat(' ', numsep)// &
            & titlestr &
            & //repeat(' ', numsep+1)//'|'
        else
            write(uf,'(*(a))') '|'//repeat(' ', numsep)// &
            & titlestr &
            & //repeat(' ', numsep)//'|'
        end if
        write(uf,'(*(a))') '+', repeat('=', width-2), '+'

        return
    end subroutine title 

    subroutine section(uf, str, sep, width)
        implicit none
        integer, intent(in) :: uf, width
        character(len=*), intent(in) :: str, sep
        
        call separator(uf, width)
        call center(uf, str, sep, width)
        call separator(uf, width)

        return
    end subroutine section 

end module strings
