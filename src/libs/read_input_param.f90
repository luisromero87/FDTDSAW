!     
! File:   read_input_param.f90
! Author: ludwig
!
! Created on April 12, 2014, 9:35 PM
!

SUBROUTINE read_input_param()
    USE Type_Kinds
    USE Constants_Module
    USE global_vars
    IMPLICIT NONE

    INTEGER :: st

    CHARACTER (LEN = name_len) :: param = ''

    OPEN(UNIT = 12, FILE = sim_config, IOSTAT = st)

    DO WHILE (param /= 'Material:')
        READ(12, *) param, material
    END DO
    
    DO WHILE (param /= 'Nstep:')
        READ(12, *) param, Nstep
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
    DO WHILE (param /= 'dt:')
        READ(12, *) param, dt
    END DO
    
    DO WHILE (param /= 'dx:')
        READ(12, *) param, deltax
    END DO
    
    DO WHILE (param /= 'dy:')
        READ(12, *) param, deltay
    END DO
    
    DO WHILE (param /= 'dz:')
        READ(12, *) param, deltaz
    END DO
    
    DO WHILE (param /= 'Nprocsx:')
        READ(12, *) param, Nprocsx
    END DO
    
    DO WHILE (param /= 'Nprocsy:')
        READ(12, *) param, Nprocsy
    END DO
    
    
    NCeldas=Nx*Ny*Nz

    CLOSE(12)
    WRITE(*,*) "\n\n"
    WRITE(*,'(A)') "**** INPUT PARAMETERS ****\n"
    WRITE(*,'(A,A)') "Material:\t", material
    WRITE(*,'(A,I7.0)') "\nNx:\t", Nx
    WRITE(*,'(A,I7.0)') "Ny:\t", Ny
    WRITE(*,'(A,I7.0)') "Nz:\t", Nz
    WRITE(*,'(A,I7.0)') "NCeldas:\t", NCeldas
    WRITE(*,'(A,e16.8)') "\ndx:\t", deltax
    WRITE(*,'(A,e16.8)') "dy:\t", deltay
    WRITE(*,'(A,e16.8)') "dz:\t", deltaz
    WRITE(*,'(A,I7.0)') "\nNprocsx:\t", Nprocsx
    WRITE(*,'(A,I7.0)') "Nprocsy:\t", Nprocsy
    WRITE(*,*) "\n\n"

END SUBROUTINE read_input_param



