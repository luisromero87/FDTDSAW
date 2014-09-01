PROGRAM pre_proc_d

    USE mpi
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE

    INTERFACE
        SUBROUTINE read_input_param()
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
        END SUBROUTINE read_input_param
        
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
    
    !Variables
    
    CHARACTER(LEN = name_len) :: outfile = 'D0.vtk' !Default
    CHARACTER(LEN = name_len) :: data_name = 'D0' !Default
    CHARACTER(LEN=name_len) :: which

    INTEGER :: st

    INTEGER :: axis, ix, iy, iz
    INTEGER :: procix, prociy
    INTEGER :: ND0x, ND0y, ND0z
    INTEGER(Long) :: thiscell
    REAL(Double), DIMENSION(3) :: junk
    REAL(Double), DIMENSION(:,:,:,:), POINTER :: read_D0
            
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
    
    !ALLOCATE THE MEMORY FOR D0
    ALLOCATE(D0x(0:NCeldas - 1))
    ALLOCATE(D0y(0:NCeldas - 1))
    ALLOCATE(D0z(0:NCeldas - 1))
    
    ! DEFAULT VALUES
    D0x = 0.0_dp
    D0y = 0.0_dp
    D0z = 0.0_dp
    ND0x = 420
    ND0y = 420
    ND0z = 32
        
    !TWO DIMENSIONAL (x,y) PROCESS INDEX AND OFFSET
    CALL ROLLPROC()
    offsetx=procsx*(Nx-3)*deltax
    offsety=procsy*(Ny-3)*deltay
    
    IF (me==0) THEN
        
        OPEN(UNIT=10,FILE=idt_file,STATUS='old',ACTION='read',IOSTAT=st)
        ALLOCATE(read_D0(0:2, 0:ND0x-1 , 0:ND0y-1 , 0:ND0z-1))
        DO iz = 0 , ND0z-1
        DO iy = 0 , ND0y-1
        DO ix = 0 , ND0x-1
            READ(10,*) junk, read_D0(0:2,ix,iy,iz)
        END Do
        END DO
        END DO
    
        CALL SYSTEM('mkdir -p outputdata/IDT')
        
        DO prociy=0, Nprocsy-1
        DO procix=0, Nprocsx-1
            DO iz=0, Nz-1
            DO iy=1, Ny-2
            DO ix=1, Nx-2
            thiscell=UROLL3(ix, iy, iz)
                D0x(thiscell)=read_D0(0,ix-1+procix*(Nx-2),iy-1+prociy*(Ny-2),iz)
                D0y(thiscell)=read_D0(1,ix-1+procix*(Nx-2),iy-1+prociy*(Ny-2),iz)
                D0z(thiscell)=read_D0(2,ix-1+procix*(Nx-2),iy-1+prociy*(Ny-2),iz)
            END DO
            END DO
            END DO
            
            offsetx=procix*(Nx-3)*deltax
            offsety=prociy*(Ny-3)*deltay
            
            WRITE(outfile, '(A,i3.3)') 'D0',procix+Nprocsx*prociy
            OPEN(UNIT=12, FILE='outputdata/IDT/'//outfile, ACTION="write", STATUS="replace", FORM = 'unformatted', IOSTAT = st)
            DO ix = 0, Nx - 1
            WRITE(12) ((D0x(UROLL3(ix, iy, iz)), iy = 0, Ny - 1), iz = 0, Nz - 1)
            WRITE(12) ((D0y(UROLL3(ix, iy, iz)), iy = 0, Ny - 1), iz = 0, Nz - 1)
            WRITE(12) ((D0z(UROLL3(ix, iy, iz)), iy = 0, Ny - 1), iz = 0, Nz - 1)
            END DO
            CLOSE(12)
            
            WRITE(outfile, '(A,i3.3,A)') 'D0_',procix+Nprocsx*prociy,'.vtk' 
            2000 format(a)
            3000 format('DIMENSIONS ',i4,i4,i4)
            4000 format(1p,e12.6)
                        
            OPEN(UNIT=12, FILE='outputdata/IDT/'//outfile, ACTION="write", STATUS="replace")
            
            WRITE(12,2000) "# vtk DataFile Version 2.0"
            WRITE(12,2000) "Volume data"
            WRITE(12,2000) "ASCII"
            WRITE(12,2000) "DATASET RECTILINEAR_GRID"
            WRITE(12,3000) Nx-2,Ny-2,Nz-1
            
            WRITE(12,*) "X_COORDINATES", Nx-2, "double"
            DO ix = 0, Nx-3
              WRITE(12,4000) ix*deltax+offsetx
            END DO
            WRITE(12,*) "Y_COORDINATES", Ny-2, "double"
            DO iy = 0, Ny-3
              WRITE(12,4000) iy*deltay+offsety
            END DO
            WRITE(12,*) "Z_COORDINATES", Nz-1, "double"
            DO iz = 1, Nz-1
              WRITE(12,4000) iz*deltaz
            END DO
            
            5000 format(1pe10.3,1x,1pe10.3,1x,1pe10.3)
            
            WRITE(12,*) "POINT_DATA", (Nx-2)*(Ny-2)*(Nz-1)
            WRITE(12,*) "SCALARS "//data_name//" double 3"
            WRITE(12,*) "LOOKUP_TABLE default"
            DO iz = 1, Nz-1
            DO iy = 1, Ny-2
            DO ix = 1, Nx-2
            thiscell=UROLL3(ix, iy, iz)            
                write(12,5000) REAL(D0x(thiscell),Single),REAL(D0y(thiscell),Single),REAL(D0z(thiscell),Single)
            END DO
            END DO 
            END DO
            
            CLOSE(12)
        END DO
        END DO

!        2000 format(a)
!        3000 format('DIMENSIONS ',i4,i4,i4)
!        4000 format(1p,e12.6)
!        
!        OPEN(UNIT=12, FILE='outputdata/IDT/'//outfile, ACTION="write", STATUS="replace")
!        
!        WRITE(12,2000) "# vtk DataFile Version 2.0"
!        WRITE(12,2000) "Volume data"
!        WRITE(12,2000) "ASCII"
!        WRITE(12,2000) "DATASET RECTILINEAR_GRID"
!        WRITE(12,3000) ND0x-2,ND0y-2,ND0z-1
!        
!        WRITE(12,*) "X_COORDINATES", ND0x-2, "double"
!        DO ix = 0, ND0x-3
!          WRITE(12,4000) ix*deltax+offsetx
!        END DO
!        WRITE(12,*) "Y_COORDINATES", ND0y-2, "double"
!        DO iy = 0, ND0y-3
!          WRITE(12,4000) iy*deltay+offsety
!        END DO
!        WRITE(12,*) "Z_COORDINATES", ND0z-1, "double"
!        DO iz = 1, ND0z-1
!          WRITE(12,4000) iz*deltaz
!        END DO
!        
!        5000 format(1pe10.3,1x,1pe10.3,1x,1pe10.3)
!        
!        WRITE(12,*) "POINT_DATA", (ND0x-2)*(ND0y-2)*(ND0z-1)
!        WRITE(12,*) "SCALARS "//data_name//" double 3"
!        WRITE(12,*) "LOOKUP_TABLE default"
!        DO iz = 1, ND0z-1
!        DO iy = 1, ND0y-2
!        DO ix = 1, ND0x-2
!            write(12,5000) REAL(read_D0(0:2,ix,iy,iz),Single)
!        END DO
!        END DO 
!        END DO
!        
!        CLOSE(12)
    
    END IF
    
    
!    WRITE(which, '(A,I3.3)') 'outputdata/weights/weights', me
!    OPEN(UNIT = 11, FILE = which, ACTION="write", STATUS="replace", FORM = 'unformatted', IOSTAT = st)
!
!    DO axis = 1, 3
!        DO ix = 0, Nx - 1
!            WRITE(11) ((w1(UROLL3(ix, iy, iz), axis), iy = 0, Ny - 1), iz = 0, Nz - 1)
!        END DO
!    END DO
!
!    DO axis = 1, 3
!        DO ix = 0, Nx - 1
!            WRITE(11) ((w2(UROLL3(ix, iy, iz), axis), iy = 0, Ny - 1), iz = 0, Nz - 1)
!        END DO
!    END DO
!
!    CLOSE(11)
!	
!	
!    CALL load_PML() 
!    WRITE(outfile, '(A,i3.3,A)') 'w1_',me,'.vtk' 
!    CALL open_vtk_file(outfile)
!    CALL write_volume_w1(outfile, data_name)
    
        
    !DEALLOCATE ALL VARIABLES
    DEALLOCATE(D0z)
    DEALLOCATE(D0y)
    DEALLOCATE(D0x)
    IF (me .EQ. 0) THEN
        DEALLOCATE(read_D0)
    END IF
    ! Close out MPI
    CALL mpi_finalize(ierr)



END PROGRAM
