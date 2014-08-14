!     
! File:   set_up_sim.f90
! Author: ludwig
!
! Created on April 12, 2014, 9:35 PM
!

SUBROUTINE set_up_sim(material, NCeldas,Nx,Ny,Nz)
    USE Type_Kinds
    USE Constants_Module
    IMPLICIT NONE
    
    CHARACTER(LEN=name_len) :: material
    INTEGER(Long) :: NCeldas
    INTEGER(Long) :: Nx
    INTEGER(Long) :: Ny
    INTEGER(Long) :: Nz

    INTEGER :: st

    CHARACTER (LEN = name_len) :: param = ''

    OPEN(UNIT = 12, FILE = sim_config, IOSTAT = st)

    DO WHILE (param /= 'Material:')
        READ(12, *) param, material
    END DO
    
    DO WHILE (param /= 'Nx:')
        READ(12, *) param, Nx
    END DO
    
    DO WHILE (param /= 'Ny:')
        READ(12, *) param, Ny
    END DO
    
    DO WHILE (param /= 'Nz:')
        READ(12, *) param, Nz
    END DO
    
    NCeldas=Nx*Ny*Nz

    CLOSE(12)

    write(*,*) Nx, Ny, Nz,NCeldas

END SUBROUTINE set_up_sim



