PROGRAM acousticwaves

    USE mpi
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    USE Lib_FDTD_SAW
    USE time_step
    IMPLICIT NONE

    !Command line argument variables
    INTEGER::narg,cptArg
    LOGICAL::lookForDebug = .FALSE.
    LOGICAL::lookForConfigfile = .FALSE.
    LOGICAL :: fileExist
    LOGICAL :: error = .FALSE.
    CHARACTER(LEN = name_len) :: NAME
    
    CHARACTER(LEN = name_len) :: outfile = 'prueba.vtk' !Default
    CHARACTER(LEN = name_len) :: data_name = 'w1' !Default
    
    CHARACTER(LEN = name_len) :: zper='False'

    !Error handle
    INTEGER :: st

    !Varibles for idexing
    INTEGER :: axis, ix, iy, iz, cont
    
    REAL(Double) :: TIME1,TIME2
       
    ! Initialize MPI; get total number of tasks and ID of this task
    CALL mpi_init(ierr)
    CALL mpi_comm_size(MPI_COMM_WORLD, ntasks, ierr)
    CALL mpi_comm_rank(MPI_COMM_WORLD, me, ierr)
    
        
    !PROCESS 0 READS THE INPUT PARAMETERS AND BROADCASTS THE INFORMATION
    IF (me==0) THEN
    !Check if arguments are found
	narg=command_argument_count()
    
    IF(narg>0)THEN
	    !loop across options
	    DO cptArg=1,narg
	        CALL get_command_argument(cptArg,NAME)
		    SELECT CASE(ADJUSTL(NAME))
			    !First known args
			    CASE("--debug")
                    Debug = "True"
                    lookForDebug = .TRUE. !change logical value
			    CASE("--configfile")
                    lookForDebug = .FALSE.
                    lookForConfigfile = .TRUE. !change logical value
			    CASE("--help")
                    WRITE(*,*) "\nAvailable options:\n\n&
                    --debug\t\t Print adicional information about the simulation\n&
                    --configfile\t input configuration file\n&
                    --help\t\t print this message\n"
                    STOP
	            CASE DEFAULT
	            !Treat the second arg of a serie
	            IF (LookForDebug) THEN
			        Debug = ADJUSTL(NAME) !assign a value to Debug
			        LookForDebug = .FALSE. !put the logical variable to its initial value
                ELSEIF (LookForConfigfile) THEN
                    input_param_file =  ADJUSTL(NAME)
                    INQUIRE (FILE = input_param_file, EXIST = fileExist)
                    IF ( .NOT. fileExist ) THEN
                        WRITE(*,*)'file ',input_param_file,' not found'
                        STOP
                    ENDIF
                    LookForConfigfile = .FALSE.
			    ELSE
			        WRITE(*,*)"Ignoring unknown option ",ADJUSTL(NAME)
			    ENDIF
		    ENDSELECT 
	    ENDDO 
	    ENDIF
    ENDIF
    CALL mpi_bcast(Debug, name_len, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    IF (me .EQ. 0) CALL read_input_param(input_param_file) 
    
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
            WRITE(*,*) "***Error: The requiered processes are ", Nprocsx*Nprocsy
        END IF
        CALL mpi_finalize(ierr)
        STOP
    END IF
    
    CALL SYSTEM('mkdir -p outputdata/IDT')
    
    !ALLOCATE THE MEMORY FOR ALL VARIABLES
    CALL allocate_memory()
    
    !PROCESS 0 READS MATERIAL CONSTANTS AND BROADCASTS THE INFORMATION
    IF (me==0) THEN
        CALL load_material(Debug=Debug)    !ok
    END IF
    CALL mpi_bcast(rho_inv, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(c_E, 36, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(s_E, 36, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(beta_s, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(e_piezo, 18, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(PMLwidth, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(smax, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL mpi_bcast(M, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    rho=1.0_dp/rho_inv
!    beta_s=0.0_dp; e_piezo=0.0_dp

!    CALL read_idt('outputdata/IDT/idt2.txt',420,420,32)
    
    !SIZE OF MPI BUFFERS
    vbuffsizex=(3*Nz*Ny)
    vbuffsizey=(3*Nz*Nx)
    
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    CALL SETUP_MPI_VARS(Debug=Debug)
    CALL allocate_memory_pml()
    CALL PML_weights() !w1, w2 ...ok
!    CALL load_D0() !w1, w2 ...ok
    WRITE(outfile, '(A,i3.3,A)') 'IDT/D0',me,'.vtr'
!    CALL save_D0_to_vtk(outfile, nx1=1, nx2=Nx-2 , ny1=1, ny2=Ny-2, nz1=1, nz2=Nz-2)
!    CALL write_volume_D0()
!    IF (Debug .EQ. 'True') THEN
!        CALL write_volume_w1() !ok
!        CALL write_volume_D0() !ok
!    ENDIF
    
!     CALL CPU_TIME(TIME1) 
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    TIME1=MPI_WTIME()
    
    OPEN(UNIT=31, FILE='outputdata/U_total.dat', ACTION="write", STATUS="replace")

    DO step = 1, Nstep
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL share_T()
        CALL v_half_step(zper=zper)
        CALL dot_source()
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL share_v()
        CALL T_half_step(zper=zper)
        CALL Get_Total_Kinetic_Energy()
        CALL Get_Total_Strain_Energy()
        CALL Get_Total_Electric_Energy()
        IF (MOD(STEP, 2) .EQ. 0) THEN
            1000 format(i3.3, '_'i3.3, '.vtr')
            WRITE(outfile, 1000) me, step/2
            CALL free_surface_vtk('free_surface'//outfile)
            CALL EA_to_vtk('xzplane'//outfile, nx1=1, nx2=Nx-2 , ny1=1, ny2=Ny-2, nz1=1, nz2=Nz-2)
        END IF
        IF (me==0) THEN
            WRITE (*,'(A,I5.5,A,I5.5)',advance="no") "\r ", step, "  out of  ",Nstep
            write(31,*) U_k_total, U_s_total, U_e_total
        END IF
    END DO
    
    close(31)
    
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
