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
    
    CHARACTER(LEN=name_len) :: which

    INTEGER(Long) UROLL3

    INTEGER :: st

    INTEGER :: ix, iy, iz, axis
    
    WRITE(which, '(A,I3.3)') 'outputdata/weights/weights', me
    OPEN(UNIT = 11, FILE = which, FORM = 'unformatted', IOSTAT = st)

    DO axis = 1, 3
        !write(*, *) axis
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

SUBROUTINE load_D02()
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE
    
    CHARACTER(LEN=name_len) :: which

    INTEGER(Long) UROLL3

    INTEGER :: st

    INTEGER :: ix, iy, iz, axis
    
    WRITE(which, '(A,I3.3)') 'outputdata/IDT/D0', me
    OPEN(UNIT = 12, FILE = which, ACTION="read", STATUS="old", FORM = 'unformatted', IOSTAT = st)
    
    IF (me .EQ. 0 .AND. st .EQ. 0) THEN
        DO ix = 0, Nx - 1
        READ(12) ((D0x(UROLL3(ix, iy, iz)), iy = 0, Ny - 1), iz = 0, Nz - 1)
        READ(12) ((D0y(UROLL3(ix, iy, iz)), iy = 0, Ny - 1), iz = 0, Nz - 1)
        READ(12) ((D0z(UROLL3(ix, iy, iz)), iy = 0, Ny - 1), iz = 0, Nz - 1)
        END DO
    ELSE IF (me .EQ. 0) THEN
        WRITE(*,*) "No input data for D0, setting D0 to 0.0"
    END IF
    
    CLOSE(12)

END SUBROUTINE load_D02

