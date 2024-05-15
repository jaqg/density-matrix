! +--------------------------------------------+
! | Author: Jose Antonio Quinonero Gris        |
! | Creation date: Tuesday 17:27:38 30/04/2024 |
! +--------------------------------------------+
program main
    !
    ! Program to compute the HF energy from the 1e and 2e-density matrices
    !
    ! Modules
    use declarations
    use strings
    use IO
    use density_matrices
    !
    implicit none
    !
    ! === START OF THE PROGRAM ===
    !
     
    ! Initialize the parameters
    call parameters

    ! Open files
    ! call open_files

    ! Print the title
    call title(6, 'Density Matrices', termwidth)
    ! call title(ouf, 'Program AOs', termwidth)


    !
    stop
endprogram main
