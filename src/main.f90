PROGRAM acousticwaves

    USE mpi
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE

    INTERFACE

        SUBROUTINE load_material()
            USE Type_Kinds
            USE Constants_Module
            IMPLICIT NONE
            CHARACTER(LEN = name_len) :: material
            REAL(Double) :: rho_inv
            REAL(Double), DIMENSION(:), POINTER :: c_E
            REAL(Double), DIMENSION(:), POINTER :: beta_s
            REAL(Double), DIMENSION(:), POINTER :: e_piezo
        END SUBROUTINE load_material

        SUBROUTINE ROLLPROC()
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
        END SUBROUTINE ROLLPROC

        SUBROUTINE load_PML()
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
        END SUBROUTINE load_PML

        SUBROUTINE load_D0()
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
        END SUBROUTINE load_D0
        
        SUBROUTINE share_v()
	        USE MPI
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
        END SUBROUTINE share_v
        
        SUBROUTINE share_T()
	        USE MPI
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
        END SUBROUTINE share_T

        SUBROUTINE read_input_param()
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
        END SUBROUTINE read_input_param

        SUBROUTINE allocate_memory()
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
        END SUBROUTINE allocate_memory

        SUBROUTINE deallocate_memory()
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
        END SUBROUTINE deallocate_memory

        SUBROUTINE open_vtk_file(outfile)
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE

            CHARACTER(LEN = name_len), INTENT(IN) :: outfile
        END SUBROUTINE open_vtk_file

        SUBROUTINE write_free_surface(outfile, data_name)
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
            CHARACTER(LEN = name_len), INTENT(IN) :: outfile
            CHARACTER(LEN = name_len), INTENT(IN) :: data_name
        END SUBROUTINE write_free_surface

        SUBROUTINE write_volume_v(outfile, data_name)
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
            CHARACTER(LEN = name_len), INTENT(IN) :: outfile
            CHARACTER(LEN = name_len), INTENT(IN) :: data_name
        END SUBROUTINE write_volume_v

        SUBROUTINE write_volume_w1()
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
        END SUBROUTINE write_volume_w1

        SUBROUTINE write_volume_D0()
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE 
        END SUBROUTINE write_volume_D0

        SUBROUTINE xzplane()
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE 
        END SUBROUTINE xzplane
    END INTERFACE

    !Functions
    INTEGER(Long) :: UROLL3
    INTEGER(Short) :: UROLLPROC
!     DOUBLE PRECISION MPI_WTIME()
    
    CHARACTER(LEN = name_len) :: outfile = 'prueba.vtk' !Default
    CHARACTER(LEN = name_len) :: data_name = 'w1' !Default

    INTEGER :: st

    INTEGER :: axis, ix, iy, iz, cont
    
    REAL(Double) :: TIME1,TIME2
    
    ! Initialize MPI; get total number of tasks and ID of this task
    CALL mpi_init(ierr)
    CALL mpi_comm_size(MPI_COMM_WORLD, ntasks, ierr)
    CALL mpi_comm_rank(MPI_COMM_WORLD, me, ierr)
    
    !PROCESS 0 READS THE INPUT PARAMETERS AND BROADCASTS THE INFORMATION
    IF (me==0) THEN
        CALL read_input_param() 
    END IF
    CALL mpi_bcast(material, name_len, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(Nstep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(Nx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(Ny, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(Nz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(NCeldas, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(deltax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(deltay, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(deltaz, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(Nprocsx, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(Nprocsy, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    IF (Nprocsx*Nprocsy .NE. ntasks) THEN
        IF (me==0) THEN
            WRITE(*,*) "The requiered process is ", Nprocsx*Nprocsy
        END IF
        CALL mpi_finalize(ierr)
        STOP
    END IF
    
    CALL SYSTEM('mkdir -p outputdata/IDT')
    
    !ALLOCATE THE MEMORY FOR ALL VARIABLES
    CALL allocate_memory()
    
    !PROCESS 0 READS MATERIAL CONSTANTS AND BROADCASTS THE INFORMATION
    IF (me==0) THEN
        CALL load_material()    !ok
    END IF
    CALL mpi_bcast(rho_inv, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(c_E, 36, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(beta_s, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(e_piezo, 18, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    
    !SIZE OF MPI BUFFERS
    vbuffsizex=(3*Nz*Ny)
    vbuffsizey=(3*Nz*Nx)
    
    !TWO DIMENSIONAL (x,y) PROCESS INDEX AND OFFSET
    CALL ROLLPROC()
    offsetx=procsx*(Nx-3)*deltax
    offsety=procsy*(Ny-3)*deltay
    !TWO DIMENSIONAL FIRST NEIGHBORS
    nUP = INT(UROLLPROC(procsx, procsy+1))
    nDOWN = INT(UROLLPROC(procsx, procsy-1))
    nRIGHT = INT(UROLLPROC(procsx+1, procsy))
    nLEFT = INT(UROLLPROC(procsx-1, procsy))
    
    !Testing ROLLPROC and UROLLPROC
!     DO ix=0, ntasks
!         IF (me==ix) THEN
!             write(*,*) procsx, procsy, UROLLPROC(procsx, procsy), nUP, nDOWN,nRIGHT,nLEFT
!         END IF
!         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!     END DO
    
    CALL load_PML() !w1, w2 ...ok
    CALL load_D0() !w1, w2 ...ok
    CALL write_volume_w1() !ok
    CALL write_volume_D0() !ok
    
!     CALL CPU_TIME(TIME1) 
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    TIME1=MPI_WTIME()
    


    DO step = 1, Nstep
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL share_T()
        CALL v_half_step()
        CALL free_boundary_v()
!        CALL dot_source()
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL share_v()
        CALL T_half_step()
        CALL free_boundary_T()
        IF (MOD(STEP, 100) .EQ. 0) THEN
            CALL xzplane()
            1000 format('free_surface', i3.3, '_'i3.3, '.vtk')
            WRITE(outfile, 1000) me, step/100
            data_name = 'v'
            CALL write_free_surface(outfile, data_name)
        END IF
!         cont = 1
!         DO ix = Nx/2-16, Nx/2+16, 4
!             DO iy = Ny/2-16, Ny/2+16, 4
!                 DO iz = 1, 2*8+1, 2
!                     probe(:, step, cont) = (/ Vx(UROLL3(ix, iy, iz)), Vy(UROLL3(ix, iy, iz)), Vz(UROLL3(ix, iy, iz)) /)
!                     cont = cont + 1
!                 END DO
!             END DO
!         END DO
        IF (me==0) THEN
            WRITE (*,'(A,I5.5,A,I5.5)',advance="no") "\r ", step, "  out of  ",Nstep
        END IF
    END DO
    
!     CALL CPU_TIME(TIME2) 
    TIME2 = MPI_WTIME()
    TIME1 = TIME2-TIME1
    
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    CALL MPI_REDUCE(TIME1,TIME2,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    IF (me .EQ. 0) THEN
        WRITE(*,*) TIME2
        OPEN(UNIT=12, FILE='time', ACTION="write", position="append")
        WRITE(12,*) TIME2
    END IF
    
    CALL deallocate_memory()
    ! Close out MPI
    CALL mpi_finalize(ierr)




END PROGRAM
