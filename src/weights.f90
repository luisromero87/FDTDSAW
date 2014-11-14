REAL(Double) FUNCTION fomega(ix)
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    USE Lib_FDTD_SAW
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ix
    
    fomega = smax*((ix+1.0_dp)/(PMLwidth*1.0_dp))**m

    RETURN

END FUNCTION fomega

PROGRAM weights

    USE mpi
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    USE Lib_FDTD_SAW
    IMPLICIT NONE

    INTERFACE
!        SUBROUTINE read_input_param()
!            USE Type_Kinds
!            USE Constants_Module
!            USE Global_Vars
!            IMPLICIT NONE
!        END SUBROUTINE read_input_param
        
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

        SUBROUTINE open_vtk_file(outfile)
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE

            CHARACTER(LEN = name_len), INTENT(IN) :: outfile
        END SUBROUTINE

        SUBROUTINE write_volume_w1()
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
        END SUBROUTINE 
        
        REAL(Double) FUNCTION fomega(ix)
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            USE Lib_FDTD_SAW
            IMPLICIT NONE

            INTEGER, INTENT(IN) :: ix
        END FUNCTION fomega
    END INTERFACE

    !Functions
    INTEGER(Long) :: UROLL3
!    REAL(Double) :: fomega
    
    !Variables
    
    CHARACTER(LEN = name_len) :: outfile = 'prueba.vtk' !Default
    CHARACTER(LEN = name_len) :: data_name = 'w1' !Default
    CHARACTER(LEN=name_len) :: which

    INTEGER :: st

    INTEGER :: axis, ix, iy, iz
    
    
    REAL(Double) :: omega
        
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
            WRITE(*,*) "Min procs requiered is ", Nprocsx*Nprocsy
        END IF
        CALL mpi_finalize(ierr)
        STOP
    END IF
    
    !ALLOCATE THE MEMORY FOR ALL VARIABLES
    ALLOCATE(w1(0:NCeldas - 1, 1:3))
    ALLOCATE(w2(0:NCeldas - 1, 1:3))
    
    ! DEFAULT VALUES
    w1 = 1.0_dp  
    w2 = dt/deltax
        
    m = 1.55_dp
    PMLwidth = 0
    smax=0.0041_dp*PMLwidth/dt
        
    !TWO DIMENSIONAL (x,y) PROCESS INDEX AND OFFSET
    CALL ROLLPROC()
    offsetx=procsx*(Nx-3)*deltax
    offsety=procsy*(Ny-3)*deltay
    
    DO ix=0,PMLwidth-1
		omega=fomega(ix)
	    DO iz=0,Nz-1
	    DO iy=0,Ny-1
			!<- x
			IF (procsx==0) THEN
				w1(UROLL3(PMLwidth-1-ix,iy,iz),1) = (2.0_dp-dt*omega)/(2.0_dp+dt*omega)
				w2(UROLL3(PMLwidth-1-ix,iy,iz),1) = 2.0_dp*dt /( deltax*(2.0_dp+dt*omega) )
			END IF
			!   x ->
			IF (procsx==Nprocsx-1) THEN
				w1(UROLL3(Nx-PMLwidth+ix,iy,iz),1) = (2.0_dp-dt*omega)/(2.0_dp+dt*omega)
				w2(UROLL3(Nx-PMLwidth+ix,iy,iz),1) = 2.0_dp*dt /( deltax*(2.0_dp+dt*omega) )
			END IF
	    END DO
	    END DO
    END DO
    
    DO iy=0,PMLwidth-1
		omega=fomega(iy)
	    DO iz=0,Nz-1
	    DO ix=0,Nx-1
			!<- y
			IF (procsy .EQ. 0) THEN
				w1(UROLL3(ix,PMLwidth-1-iy,iz),2) = (2.0_dp-dt*omega)/(2.0_dp+dt*omega)
				w2(UROLL3(ix,PMLwidth-1-iy,iz),2) = 2.0_dp*dt /( deltax*(2.0_dp+dt*omega) )
			END IF
			!   y ->
			IF (procsy .EQ. Nprocsy-1) THEN
				w1(UROLL3(ix,Ny-PMLwidth+iy,iz),2) = (2.0_dp-dt*omega)/(2.0_dp+dt*omega)
				w2(UROLL3(ix,Ny-PMLwidth+iy,iz),2) = 2.0_dp*dt /( deltax*(2.0_dp+dt*omega) )
			END IF
	    END DO
	    END DO
    END DO
    
    DO iz=0,PMLwidth-1
		omega=fomega(iz)
	    DO iy=0,Ny-1
	    DO ix=0,Nx-1
			!   z ->
				w1(UROLL3(ix,iy,Nz-PMLwidth+iz),3) = (2.0_dp-dt*omega)/(2.0_dp+dt*omega)
				w2(UROLL3(ix,iy,Nz-PMLwidth+iz),3) = 2.0_dp*dt /( deltax*(2.0_dp+dt*omega) )
	    END DO
	    END DO
    END DO
    

	CALL SYSTEM('mkdir -p outputdata/weights')
    
    WRITE(which, '(A,I3.3)') 'outputdata/weights/weights', me
    OPEN(UNIT = 11, FILE = which, ACTION="write", STATUS="replace", FORM = 'unformatted', IOSTAT = st)

    DO axis = 1, 3
        !write(*, *) axis
        DO ix = 0, Nx - 1
            WRITE(11) ((w1(UROLL3(ix, iy, iz), axis), iy = 0, Ny - 1), iz = 0, Nz - 1)
        END DO
    END DO

    DO axis = 1, 3
        DO ix = 0, Nx - 1
            WRITE(11) ((w2(UROLL3(ix, iy, iz), axis), iy = 0, Ny - 1), iz = 0, Nz - 1)
        END DO
    END DO

    CLOSE(11)
	
	
    CALL load_PML() 
    CALL write_volume_w1()
    
        
    !DEALLOCATE ALL VARIABLES
    DEALLOCATE(w1)
    DEALLOCATE(w2)
    ! Close out MPI
    CALL mpi_finalize(ierr)



END PROGRAM
