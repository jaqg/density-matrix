! +--------------------------------------------+
! | Author: Jose Antonio Quinonero Gris        |
! | Creation date: Tuesday 09:11:31 16/01/2024 |
! +--------------------------------------------+
module declarations 
    ! 
    ! Module to declare variables
    !
    implicit none
    integer :: termwidth 

    contains

    subroutine parameters
        implicit none
        
        ! Width of the terminal/file to be written on
        termwidth = 46

        return
    end subroutine parameters 

    !
end module declarations 
