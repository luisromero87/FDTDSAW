INTEGER(Long) FUNCTION UROLL3(ix, iy, iz)
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE

    INTEGER(Long), INTENT(IN) :: ix, iy, iz

    UROLL3 = MOD(iz + Nz, Nz) * Ny * Nx + MOD(iy + Ny, Ny) * Nx + MOD(ix + Nx, Nx)

    RETURN

END FUNCTION UROLL3

