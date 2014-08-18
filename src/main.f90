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
        END SUBROUTINE

        SUBROUTINE write_volume_v(outfile, data_name)
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
            CHARACTER(LEN = name_len), INTENT(IN) :: outfile
            CHARACTER(LEN = name_len), INTENT(IN) :: data_name
        END SUBROUTINE

        SUBROUTINE write_volume_w1(outfile, data_name)
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
            CHARACTER(LEN = name_len), INTENT(IN) :: outfile
            CHARACTER(LEN = name_len), INTENT(IN) :: data_name
        END SUBROUTINE 
    END INTERFACE

    !Functions
    INTEGER(Long) :: UROLL3
    INTEGER(Short) :: UROLLPROC
    
    CHARACTER(LEN = name_len) :: outfile = 'prueba.vtk' !Default
    CHARACTER(LEN = name_len) :: data_name = 'w1' !Default

    INTEGER :: st

    INTEGER :: axis, ix, iy, iz, cont

    REAL(Double), DIMENSION(:,:,:), POINTER :: probe
    !REAL(Double), DIMENSION(1:3,1:Nstep, 1:Nprobe) :: probe = 0.0_dp
    ALLOCATE(probe(1:3,1:Nstep, 1:Nprobe))
    probe=0.0_dp
    
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
    IF (Nprocsx*Nprocsy>ntasks) THEN
        IF (me==0) THEN
            WRITE(*,*) "Min procs requiered is ", Nprocsx*Nprocsy
        END IF
        CALL mpi_finalize(ierr)
        STOP
    END IF
    
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
    
    !TWO DIMENSIONAL (x,y) PROCESS INDEX AND OFFSET
    CALL ROLLPROC()
    offsetx=procsx*(Nx-3)*deltax
    offsety=procsy*(Ny-3)*deltay
    
    !Testing ROLLPROC and UROLLPROC
    !DO ix=0, ntasks
    !    IF (me==ix) THEN
    !        write(*,*) procsx, procsy, UROLLPROC(procsx, procsy)
    !    END IF
    !    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !END DO
    
    CALL load_PML() !w1, w2 ...ok
    WRITE(outfile, '(A,i3.3,A)') 'w1_',me,'.vtk' 
    CALL open_vtk_file(outfile)
    CALL write_volume_w1(outfile, data_name) !ok

    DO step = 1, Nstep
        CALL v_half_step()
        CALL free_boundary_v()
        IF (me==1) THEN
	        CALL dot_source()
        END IF
        CALL share_v()
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL T_half_step()
        CALL free_boundary_T()
        CALL share_T()
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        IF (MOD(STEP, 100) .EQ. 0) THEN
            1000 format('onda', i3.3, '_'i3.3, '.vtk')
            WRITE(outfile, 1000) me, step/100
            data_name = 'v'
            CALL open_vtk_file(outfile)
            CALL write_volume_v(outfile, data_name)
        END IF
        cont = 1
        DO ix = Nx/2-16, Nx/2+16, 4
            DO iy = Ny/2-16, Ny/2+16, 4
                DO iz = 1, 2*8+1, 2
                    probe(:, step, cont) = (/ Vx(UROLL3(ix, iy, iz)), Vy(UROLL3(ix, iy, iz)), Vz(UROLL3(ix, iy, iz)) /)
                    cont = cont + 1
                END DO
            END DO
        END DO
        IF (me==0) THEN
            WRITE (*,'(A,I5.5,A,I5.5)',advance="no") "\r ", step, "  out of  ",Nstep
        END IF
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    END DO
    
    
    CALL deallocate_memory()
    ! Close out MPI
    CALL mpi_finalize(ierr)
!~     OPEN(UNIT = 12, FILE = 'probe', ACTION = "write", STATUS = "replace", form = 'unformatted')
!~     write(12) probe
!~     CLOSE(12)




END PROGRAM
