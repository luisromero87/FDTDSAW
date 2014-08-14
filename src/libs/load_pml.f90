!     
! File:   load_pml.f90
! Author: ludwig
!
! Created on April 12, 2014, 9:35 PM
!

SUBROUTINE load_PML()
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE

    INTEGER(Long) UROLL3

    INTEGER :: st

    INTEGER :: ix, iy, iz, axis

    OPEN(UNIT = 11, FILE = PLM_archive, FORM = 'unformatted', IOSTAT = st)

    DO axis = 1, 3
        write(*, *) axis
        DO ix = 0, Nx - 1
            READ(11) ((w1(UROLL3(ix, iy, iz), axis), iy = 0, Ny - 1), iz = 0, Nz - 1)
        END DO
    END DO

    DO axis = 1, 3
        DO ix = 0, Nx - 1
            READ(11) ((w2(UROLL3(ix, iy, iz), axis), iy = 0, Ny - 1), iz = 0, Nz - 1)
        END DO
    END DO

    CLOSE(11)

END SUBROUTINE load_PML

