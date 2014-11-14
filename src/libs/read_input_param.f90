!     
! File:   read_input_param.f90
! Author: ludwig
!
! Created on April 12, 2014, 9:35 PM
!

SUBROUTINE read_input_param2()
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
    
    DO WHILE (param /= 'Nstep:')
        READ(12, *) param, Nstep
    END DO
    DO WHILE (param /= 'Nx:')
        READ(12, *) param, NGx
    END DO
    DO WHILE (param /= 'Ny:')
        READ(12, *) param, NGy
    END DO
    DO WHILE (param /= 'Nz:')
        READ(12, *) param, NGz
    END DO
    
    DO WHILE (param /= 'Nprocsx:')
        READ(12, *) param, Nprocsx
    END DO
    DO WHILE (param /= 'Nprocsy:')
        READ(12, *) param, Nprocsy
    END DO
    
    NGCeldas=NGx*NGy*NGz

    CLOSE(12)
    WRITE(*,*) "\n\n"
    WRITE(*,'(A)') "**** INPUT PARAMETERS ****\n"
    WRITE(*,'(A,A)') "Material:\t", material
    
    WRITE(*,'(A,e16.8)') "\ndt:\t", dt
    WRITE(*,'(A,e16.8)') "dx:\t", deltax
    WRITE(*,'(A,e16.8)') "dy:\t", deltay
    WRITE(*,'(A,e16.8)') "dz:\t", deltaz
    
    WRITE(*,'(A,I7.0)') "\nNstep:\t", Nstep
    WRITE(*,'(A,I7.0)') "\nNx:\t", NGx
    WRITE(*,'(A,I7.0)') "Ny:\t", NGy
    WRITE(*,'(A,I7.0)') "Nz:\t", NGz
    WRITE(*,'(A,I7.0)') "NCeldas:\t", NGCeldas
    
    WRITE(*,'(A,I7.0)') "\nNprocsx:\t", Nprocsx
    WRITE(*,'(A,I7.0)') "Nprocsy:\t", Nprocsy
    WRITE(*,'(A,I7.0)') "Total procs:\t", Nprocsx*Nprocsy
    WRITE(*,*) "\n\n"

END SUBROUTINE read_input_param2



