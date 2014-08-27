INTEGER(Long) FUNCTION UROLL3(ix, iy, iz)
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE

    INTEGER(Long), INTENT(IN) :: ix, iy, iz

    UROLL3 = MOD(iz + Nz, Nz) * Ny * Nx + MOD(iy + Ny, Ny) * Nx + MOD(ix + Nx, Nx)

    RETURN

END FUNCTION UROLL3


INTEGER(Short) FUNCTION UROLLPROC(px, py)
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE

    INTEGER(Short), INTENT(IN) :: px, py

    UROLLPROC = MOD(py+Nprocsy,Nprocsy) * Nprocsx + MOD(px+Nprocsx,Nprocsx)

    RETURN

END FUNCTION UROLLPROC

SUBROUTINE ROLLPROC()
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE

    procsy = INT(me/Nprocsx,Short)
    procsx = MOD(me,Nprocsx)

END SUBROUTINE ROLLPROC

