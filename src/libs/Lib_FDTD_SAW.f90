MODULE Lib_FDTD_SAW

USE mpi
USE Type_Kinds
USE Constants_Module
USE Global_Vars
USE LIB_VTK_IO

IMPLICIT NONE
PRIVATE

PUBLIC :: read_input_param
PUBLIC :: allocate_memory
PUBLIC :: allocate_memory_pml
PUBLIC :: deallocate_memory
PUBLIC :: load_material
PUBLIC :: SETUP_MPI_VARS
PUBLIC :: PML_weights
PUBLIC :: load_D0
PUBLIC :: freesurface_to_vtk
PUBLIC :: EA_to_vtk
PUBLIC :: read_idt
PUBLIC :: save_D0_to_vtk
PUBLIC :: Get_Total_Kinetic_Energy
PUBLIC :: Get_Total_Strain_Energy
PUBLIC :: Get_Total_Electric_Energy

CONTAINS

REAL(Double) FUNCTION fomega(ix)
    IMPLICIT NONE
    INTEGER(Long), INTENT(IN) :: ix
    fomega = smax*((ix+1.0_dp)/(PMLwidth*1.0_dp))**m
    RETURN
ENDFUNCTION fomega

INTEGER(Long) FUNCTION UROLL3(ix, iy, iz)
    IMPLICIT NONE
    INTEGER(Long), INTENT(IN) :: ix, iy, iz
    UROLL3 = MOD(iz + Nz, Nz) * Ny * Nx + MOD(iy + Ny, Ny) * Nx + MOD(ix + Nx, Nx)
    RETURN
ENDFUNCTION UROLL3

INTEGER(Short) FUNCTION UROLLPROC(px, py)
    IMPLICIT NONE
    INTEGER(Short), INTENT(IN) :: px, py
    UROLLPROC = MOD(py+Nprocsy,Nprocsy) * Nprocsx + MOD(px+Nprocsx,Nprocsx)
    RETURN
ENDFUNCTION UROLLPROC

SUBROUTINE read_input_param(input_param_file)
    IMPLICIT NONE

    CHARACTER(LEN = *) :: input_param_file
        
    INTEGER :: st
    CHARACTER (LEN = name_len) :: param = ''
    CHARACTER (LEN = name_len) :: param_name = ''
    CHARACTER (LEN = name_len) :: param_value = ''
 
    OPEN(UNIT = 12, FILE = input_param_file, IOSTAT = st)
    IF (st .EQ. 0) THEN
    READ(12, *, IOSTAT=st) param_name, param_value
    DO WHILE(st .EQ. 0)
        
        SELECT CASE(param_name)
            CASE("Material:")
                material=param_value
            CASE("dt:")
                READ(param_value,*) dt
            CASE("dx:")
                READ(param_value,*) deltax
            CASE("dy:")
                READ(param_value,*) deltay
            CASE("dz:")
                READ(param_value,*) deltaz
            CASE("Nstep:")
                READ(param_value,*) Nstep
            CASE("Data_fstep:")
                READ(param_value,*) data_fstep
            CASE("Nx:")
                READ(param_value,*) NGx
            CASE("Ny:")
                READ(param_value,*) NGy
            CASE("Nz:")
                READ(param_value,*) NGz
            CASE("Nprocsx:")
                READ(param_value,*) Nprocsx
            CASE("Nprocsy:")
                READ(param_value,*) Nprocsy
            CASE("PML_width:")
                READ(param_value,*) PMLwidth
            CASE("smax:")
                READ(param_value,*) smax
            CASE("m:")
                READ(param_value,*) m
            CASE("Output_dir:")
                output_dir=TRIM(param_value)//'/'
        ENDSELECT 
        READ(12, *, IOSTAT=st) param_name, param_value
        
    ENDDO
    CLOSE(12)
    
    WRITE(*,*) "\n"
    WRITE(*,*)  "\n*************** INPUT PARAMETERS **************\n"
    WRITE(*,'(A,A)') "Material:\t\t", material
    
    WRITE(*,'(A,'//FR4P//')') "\ndt:\t", dt
    WRITE(*,'(A,'//FR4P//')') "dx:\t", deltax
    WRITE(*,'(A,'//FR4P//')') "dy:\t", deltay
    WRITE(*,'(A,'//FR4P//')') "dz:\t", deltaz
    
    WRITE(*,'(A,'//FI4P//')') "\nNstep:\t\t", Nstep
    WRITE(*,'(A,'//FI4P//')') "data_fstep:\t", data_fstep
    WRITE(*,'(A,'//FI4P//')') "\nNx:\t\t", NGx
    WRITE(*,'(A,'//FI4P//')') "Ny:\t\t", NGy
    WRITE(*,'(A,'//FI4P//')') "Nz:\t\t", NGz
    
    WRITE(*,'(A,'//FI4P//')') "\nNprocsx:\t", Nprocsx
    WRITE(*,'(A,'//FI4P//')') "Nprocsy:\t", Nprocsy
    
    WRITE(*,'(A,'//FI4P//')') "\nPML_width:\t", PMLwidth
    WRITE(*,'(A,'//FR4P//')') "smax:\t", smax
    WRITE(*,'(A,'//FR4P//')') "m:\t", m
    
    WRITE(*,'(A,A)') "\nOutput_dir:\t", output_dir
    
    Nx=CEILING(1.0_dp*NGx/Nprocsx)+2
    NGx=(Nx-2)*Nprocsx
    Ny=CEILING(1.0_dp*NGy/Nprocsy)+2
    NGy=(Ny-2)*Nprocsy
    Nz=NGz+2
    NCeldas=Nx*Ny*Nz
    NGCeldas=NGx*NGy*NGz
    smax=smax*PMLwidth/dt
    
    WRITE(*,*) "\n**************** SIZE SETTING *****************\n"
    WRITE(*,*)   "                       proc               total"
    
    WRITE(*,'(A,2'//FI8P//')') "Nx:\t", Nx, NGx
    WRITE(*,'(A,2'//FI8P//')') "Ny:\t", Ny, NGx
    WRITE(*,'(A,2'//FI8P//')') "Nz:\t", Nz, NGz
    WRITE(*,'(A,2'//FI8P//')') "Celdas:\t", NCeldas, NGCeldas
    
    WRITE(*,'(A,I7.0)') "\nTotal procs:\t", Nprocsx*Nprocsy
    WRITE(*,*) "\n"
    
    ELSE
        WRITE(*,*) 'File ', input_param_file,' not found'
    ENDIF
    
ENDSUBROUTINE read_input_param


SUBROUTINE allocate_memory()
    IMPLICIT NONE
    
    ALLOCATE(c_E(6,6), s_E(6,6), beta_s(3,3), e_piezo(3,6))
    
    ALLOCATE(dx(0:Nx - 1), dy(0:Ny - 1), dz(0:Nz - 1))

    ALLOCATE(Vx(0:Nx-1,0:Ny-1,0:Nz-1), Vy(0:Nx-1,0:Ny-1,0:Nz-1), Vz(0:Nx-1,0:Ny-1,0:Nz-1))
    ALLOCATE(T1(0:Nx-1,0:Ny-1,0:Nz-1), T2(0:Nx-1,0:Ny-1,0:Nz-1), T3(0:Nx-1,0:Ny-1,0:Nz-1))
    ALLOCATE(T4(0:Nx-1,0:Ny-1,0:Nz-1), T5(0:Nx-1,0:Ny-1,0:Nz-1), T6(0:Nx-1,0:Ny-1,0:Nz-1))

    ALLOCATE(dEx(0:Nx-1,0:Ny-1,0:Nz-1), dEy(0:Nx-1,0:Ny-1,0:Nz-1), dEz(0:Nx-1,0:Ny-1,0:Nz-1))
    ALLOCATE(Ex(0:Nx-1,0:Ny-1,0:Nz-1), Ey(0:Nx-1,0:Ny-1,0:Nz-1), Ez(0:Nx-1,0:Ny-1,0:Nz-1))

    ALLOCATE(dDx(0:Nx-1,0:Ny-1,0:Nz-1), dDy(0:Nx-1,0:Ny-1,0:Nz-1), dDz(0:Nx-1,0:Ny-1,0:Nz-1))
    ALLOCATE(D0x(0:Nx-1,0:Ny-1,0:Nz-1), D0y(0:Nx-1,0:Ny-1,0:Nz-1), D0z(0:Nx-1,0:Ny-1,0:Nz-1))
    ALLOCATE(Disx(0:Nx-1,0:Ny-1,0:Nz-1), Disy(0:Nx-1,0:Ny-1,0:Nz-1), Disz(0:Nx-1,0:Ny-1,0:Nz-1))
    
    ALLOCATE(w1(0:Nx-1,0:Ny-1,0:Nz-1, 1:3), w2(0:Nx-1,0:Ny-1,0:Nz-1, 1:3))
    
    ALLOCATE(mpibufferx(0:6*Nz*Ny-1), mpibuffery(0:6*Nz*Nx-1))
    
    ALLOCATE(recvbuff_UP(0:6*Nz*Nx-1), recvbuff_DOWN(0:6*Nz*Nx-1))
    ALLOCATE(recvbuff_RIGHT(0:6*Nz*Ny-1), recvbuff_LEFT(0:6*Nz*Ny-1))
    ALLOCATE(sendbuff_UP(0:6*Nz*Nx-1), sendbuff_DOWN(0:6*Nz*Nx-1))
    ALLOCATE(sendbuff_RIGHT(0:6*Nz*Ny-1), sendbuff_LEFT(0:6*Nz*Ny-1))
    
    c_E=0.0_dp; s_E=0.0_dp; beta_s=0.0_dp; e_piezo=0.0_dp
    
    w1 = 1.0_dp; w2 = dt/deltax
    
    dx = deltax
    dy = deltay
    dz = deltaz
    
    Vx = 0.0_dp; Vy = 0.0_dp; Vz = 0.0_dp
    T1 = 0.0_dp; T2 = 0.0_dp; T3 = 0.0_dp; T4 = 0.0_dp; T5 = 0.0_dp; T6 = 0.0_dp

    dEx = 0.0_dp; dEy = 0.0_dp; dEz = 0.0_dp
    dDx = 0.0_dp; dDy = 0.0_dp; dDz = 0.0_dp
    D0x = 0.0_dp; D0y = 0.0_dp; D0z = 0.0_dp
    
ENDSUBROUTINE allocate_memory


SUBROUTINE allocate_memory_pml()
    IMPLICIT NONE
    
    INTEGER :: pmlxstart,  pmlxend,  pmlystart,  pmlyend,  pmlzstart,  pmlzend
    pmlxstart=0;  pmlxend=Nx-1;  pmlystart=0;  pmlyend=Ny-1;  pmlzstart=Nz-PMLwidth;  pmlzend=Nz-1
    
    xstart=1; xend=Nx-2; ystart=1; yend=Ny-2; zstart=1; zend=Nz-2;
    IF (PMLwidth .GT. 0) THEN
        IF (procsx .EQ. 0) THEN
            xstart=PMLwidth
        ENDIF
        IF (procsx .EQ. Nprocsx-1) THEN
            xend=Nx-PMLwidth-1
        ENDIF
        IF (procsy .EQ. 0) THEN
            ystart=PMLwidth
        ENDIF
        IF (procsy .EQ. Nprocsy-1) THEN
            yend=Ny-PMLwidth-1
        ENDIF
        zend=Nz-PMLwidth-1;
        
        IF (procsx .EQ. 0 .OR. procsx .EQ. Nprocsx-1 .OR. procsy .EQ. 0 .OR. procsy .EQ. Nprocsy-1) THEN
            pmlzstart=0
            pmlzend=Nz-1
        ENDIF
            
        ALLOCATE(Vx_x(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        ALLOCATE(Vx_y(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        ALLOCATE(Vx_z(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        
        ALLOCATE(Vy_x(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        ALLOCATE(Vy_y(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        ALLOCATE(Vy_z(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        
        ALLOCATE(Vz_x(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        ALLOCATE(Vz_y(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        ALLOCATE(Vz_z(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        
        ALLOCATE(T1_x(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        ALLOCATE(T1_y(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        ALLOCATE(T1_z(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        
        ALLOCATE(T2_x(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        ALLOCATE(T2_y(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        ALLOCATE(T2_z(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        
        ALLOCATE(T3_x(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        ALLOCATE(T3_y(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        ALLOCATE(T3_z(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        
        ALLOCATE(T4_y(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        ALLOCATE(T4_z(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        
        ALLOCATE(T5_x(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        ALLOCATE(T5_z(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        
        ALLOCATE(T6_x(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        ALLOCATE(T6_y(0:Nx-1, 0:Ny-1, pmlzstart:pmlzend))
        
        Vx_x = 0.0_dp; Vy_x = 0.0_dp; Vz_x = 0.0_dp;
        T1_x = 0.0_dp; T2_x = 0.0_dp; T3_x = 0.0_dp;
        T5_x = 0.0_dp; T6_x = 0.0_dp;
        
        Vx_y = 0.0_dp; Vy_y = 0.0_dp; Vz_y = 0.0_dp;
        T1_y = 0.0_dp; T2_y = 0.0_dp; T3_y = 0.0_dp;
        T4_y = 0.0_dp; T6_y = 0.0_dp;
        
        Vx_z = 0.0_dp; Vy_z = 0.0_dp; Vz_z = 0.0_dp;
        T1_z = 0.0_dp; T2_z = 0.0_dp; T3_z = 0.0_dp;
        T4_z = 0.0_dp; T5_z = 0.0_dp;
        
    ENDIF
    
ENDSUBROUTINE allocate_memory_pml

SUBROUTINE deallocate_memory()
    IMPLICIT NONE
    
    INTEGER :: st
        
    DEALLOCATE(c_E,STAT=st)
    DEALLOCATE(beta_s,STAT=st)
    DEALLOCATE(e_piezo,STAT=st)
    
    DEALLOCATE(dx,STAT=st)
    DEALLOCATE(dy,STAT=st)
    DEALLOCATE(dz,STAT=st)

    DEALLOCATE(Vx,STAT=st)
    DEALLOCATE(Vx_x,STAT=st)
    DEALLOCATE(Vx_y,STAT=st)
    DEALLOCATE(Vx_z,STAT=st)

    DEALLOCATE(Vy,STAT=st)
    DEALLOCATE(Vy_x,STAT=st)
    DEALLOCATE(Vy_y,STAT=st)
    DEALLOCATE(Vy_z,STAT=st)

    DEALLOCATE(Vz,STAT=st)
    DEALLOCATE(Vz_x,STAT=st)
    DEALLOCATE(Vz_y,STAT=st)
    DEALLOCATE(Vz_z,STAT=st)

    DEALLOCATE(T1,STAT=st)
    DEALLOCATE(T1_x,STAT=st)
    DEALLOCATE(T1_y,STAT=st)
    DEALLOCATE(T1_z,STAT=st)

    DEALLOCATE(T2,STAT=st)
    DEALLOCATE(T2_x,STAT=st)
    DEALLOCATE(T2_y,STAT=st)
    DEALLOCATE(T2_z,STAT=st)

    DEALLOCATE(T3,STAT=st)
    DEALLOCATE(T3_x,STAT=st)
    DEALLOCATE(T3_y,STAT=st)
    DEALLOCATE(T3_z,STAT=st)

    DEALLOCATE(T4,STAT=st)
    DEALLOCATE(T4_y,STAT=st)
    DEALLOCATE(T4_z,STAT=st)

    DEALLOCATE(T5,STAT=st)
    DEALLOCATE(T5_x,STAT=st)
    DEALLOCATE(T5_z,STAT=st)

    DEALLOCATE(T6,STAT=st)
    DEALLOCATE(T6_x,STAT=st)
    DEALLOCATE(T6_y,STAT=st)

    DEALLOCATE(dEx,STAT=st)
    DEALLOCATE(dEy,STAT=st)
    DEALLOCATE(dEz,STAT=st)

    DEALLOCATE(dDx,STAT=st)
    DEALLOCATE(dDy,STAT=st)
    DEALLOCATE(dDz,STAT=st)

    DEALLOCATE(D0x,STAT=st)
    DEALLOCATE(D0y,STAT=st)
    DEALLOCATE(D0z,STAT=st)

    DEALLOCATE(w1,STAT=st)
    DEALLOCATE(w2,STAT=st)

END SUBROUTINE deallocate_memory


SUBROUTINE load_material(Debug)
    IMPLICIT NONE
    
    CHARACTER(LEN=*), OPTIONAL :: Debug
    
    INTEGER :: st
            
    OPEN(UNIT = 10, FILE =  Materials_dir//material, FORM='unformatted',IOSTAT = st)
    
    READ(10) rho_inv
    READ(10) c_E
    READ(10) s_E
    READ(10) beta_s
    READ(10) e_piezo
    
    CLOSE(10)
    
    
    IF (PRESENT(Debug)) THEN
!    write(*,*) Debug
        IF (Debug .EQ. 'True') THEN
            WRITE(*,*) "\n************ MATERIAL CONSTANTS ***************\n"
            write(*,'(A,'//FR4P//')') 'Densidad:\t', 1.0_dp/rho_inv
            write(*,'(A)') '\nStiffnes Tensor:\n'
            write(*,'(6'//FR4P//')') c_E(1,:)
            write(*,'(6'//FR4P//')') c_E(2,:)
            write(*,'(6'//FR4P//')') c_E(3,:)
            write(*,'(6'//FR4P//')') c_E(4,:)
            write(*,'(6'//FR4P//')') c_E(5,:)
            write(*,'(6'//FR4P//')') c_E(6,:)
            write(*,'(A)') '\nPiezoelectric constants:\n'
            write(*,'(6'//FR4P//')') e_piezo(1,:)
            write(*,'(6'//FR4P//')') e_piezo(2,:)
            write(*,'(6'//FR4P//')') e_piezo(3,:)
            write(*,'(A)') '\nbeta_s:\n'
            write(*,'(3'//FR4P//')') beta_s(1,:)
            write(*,'(3'//FR4P//')') beta_s(2,:)
            write(*,'(3'//FR4P//')') beta_s(3,:)
            write(*,'(A)') '\n'
        ENDIF
    ENDIF
    
ENDSUBROUTINE load_material

SUBROUTINE SETUP_MPI_VARS(Debug)
    IMPLICIT NONE
    CHARACTER(LEN=*), OPTIONAL :: Debug
    
    INTEGER :: ix
    
    !TWO DIMENSIONAL (x,y) PROCESS INDEX AND OFFSET
    procsy = INT(me/Nprocsx,Short)
    procsx = MOD(me,Nprocsx)
    offsetx=procsx*(Nx-3)*deltax
    offsety=procsy*(Ny-3)*deltay
    !TWO DIMENSIONAL FIRST NEIGHBORS
    nUP = INT(UROLLPROC(procsx, procsy+1_Short))
    nDOWN = INT(UROLLPROC(procsx, procsy-1_Short))
    nRIGHT = INT(UROLLPROC(procsx+1_Short, procsy))
    nLEFT = INT(UROLLPROC(procsx-1_Short, procsy))
    
    IF (PRESENT(Debug)) THEN
        IF (Debug .EQ. 'True') THEN
            !Testing ROLLPROC and UROLLPROC
             DO ix=0, ntasks
                 IF (me==ix) THEN
                     write(*,*) procsx, procsy, UROLLPROC(procsx, procsy), nUP, nDOWN,nRIGHT,nLEFT
                 END IF
                 CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
             END DO
             IF (me .EQ. 0) THEN; WRITE(*,*) '\n'; ENDIF
        ENDIF
    ENDIF
    
ENDSUBROUTINE SETUP_MPI_VARS

SUBROUTINE PML_weights(Debug)
    IMPLICIT NONE
    CHARACTER(LEN=*), OPTIONAL :: Debug
    
    INTEGER :: st
    
    INTEGER(Long) :: axis, ix, iy, iz
    REAL(Double) :: omega
    
    CHARACTER(LEN = name_len) :: outfile = 'prueba.vtk' !Default
    CHARACTER(LEN = name_len) :: data_name = 'w1' !Default
    CHARACTER(LEN=name_len) :: which
    
    DO ix=0,PMLwidth-1
		omega=fomega(ix)
	    DO iz=0,Nz-1
	    DO iy=0,Ny-1
			!<- x
			IF (procsx==0) THEN
				w1(PMLwidth-1-ix,iy,iz,1) = (2.0_dp-dt*omega)/(2.0_dp+dt*omega)
				w2(PMLwidth-1-ix,iy,iz,1) = 2.0_dp*dt /( deltax*(2.0_dp+dt*omega) )
			END IF
			!   x ->
			IF (procsx==Nprocsx-1) THEN
				w1(Nx-PMLwidth+ix,iy,iz,1) = (2.0_dp-dt*omega)/(2.0_dp+dt*omega)
				w2(Nx-PMLwidth+ix,iy,iz,1) = 2.0_dp*dt /( deltax*(2.0_dp+dt*omega) )
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
				w1(ix,PMLwidth-1-iy,iz,2) = (2.0_dp-dt*omega)/(2.0_dp+dt*omega)
				w2(ix,PMLwidth-1-iy,iz,2) = 2.0_dp*dt /( deltax*(2.0_dp+dt*omega) )
			END IF
			!   y ->
			IF (procsy .EQ. Nprocsy-1) THEN
				w1(ix,Ny-PMLwidth+iy,iz,2) = (2.0_dp-dt*omega)/(2.0_dp+dt*omega)
				w2(ix,Ny-PMLwidth+iy,iz,2) = 2.0_dp*dt /( deltax*(2.0_dp+dt*omega) )
			END IF
	    END DO
	    END DO
    END DO
    
    DO iz=0,PMLwidth-1
		omega=fomega(iz)
	    DO iy=0,Ny-1
	    DO ix=0,Nx-1
			!   z ->
				w1(ix,iy,Nz-PMLwidth+iz,3) = (2.0_dp-dt*omega)/(2.0_dp+dt*omega)
				w2(ix,iy,Nz-PMLwidth+iz,3) = 2.0_dp*dt /( deltax*(2.0_dp+dt*omega) )
	    END DO
	    END DO
    END DO
    
        
    IF (PRESENT(Debug)) THEN
        IF (Debug .EQ. 'True') THEN
            CALL SYSTEM('mkdir -p '//trim(output_dir)//'weights')
            WRITE(which, '(A,I3.3)') trim(output_dir)//'weights/weights', me
            OPEN(UNIT = 11, FILE = which, ACTION="write", STATUS="replace", FORM = 'unformatted', IOSTAT = st)
            DO axis = 1, 3
                !write(*, *) axis
                DO ix = 0, Nx - 1
                    WRITE(11) ((w1(ix, iy, iz, axis), iy = 0, Ny - 1), iz = 0, Nz - 1)
                END DO
            END DO
            DO axis = 1, 3
                DO ix = 0, Nx - 1
                    WRITE(11) ((w2(ix, iy, iz, axis), iy = 0, Ny - 1), iz = 0, Nz - 1)
                END DO
            END DO
            CLOSE(11)
        ENDIF
    ENDIF
    
ENDSUBROUTINE PML_weights

SUBROUTINE load_D0()
    IMPLICIT NONE
    
    CHARACTER(LEN=name_len) :: which

    INTEGER(Long) UROLL3

    INTEGER :: st
    INTEGER :: ix, iy, iz, axis
    
    WRITE(which, '(A,I3.3)') trim(output_dir)//'IDT/D0', me
    OPEN(UNIT = 12, FILE = which, ACTION="read", STATUS="old", FORM = 'unformatted', IOSTAT = st)
    
    IF ( st .EQ. 0) THEN
        DO ix = 0, Nx - 1
        READ(12) ((D0x(ix, iy, iz), iy = 0, Ny - 1), iz = 0, Nz - 1)
        READ(12) ((D0y(ix, iy, iz), iy = 0, Ny - 1), iz = 0, Nz - 1)
        READ(12) ((D0z(ix, iy, iz), iy = 0, Ny - 1), iz = 0, Nz - 1)
        END DO
    ELSE IF (me .EQ. 0) THEN
        WRITE(*,*) "No input data for D0, setting D0 to 0.0\n"
    END IF
    
    CLOSE(12)

ENDSUBROUTINE load_D0

SUBROUTINE Get_Total_Kinetic_Energy()
	IMPLICIT NONE
    
	INTEGER :: ix, iy, iz
    INTEGER(Long) :: thiscell
    
    U_k = 0
    DO iz = zstart, zend
        DO iy = ystart, yend
            DO ix = xstart, xend
                thiscell=UROLL3(ix, iy, iz)
                U_k=U_k+Vx(ix, iy, iz)**2+Vy(ix, iy, iz)**2+Vz(ix, iy, iz)**2
            END DO
        END DO
    END DO
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(U_k, U_k_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    IF (me .EQ. 0) U_k_total=0.5*U_k_total*rho
    
ENDSUBROUTINE Get_Total_Kinetic_Energy

SUBROUTINE Get_Total_Strain_Energy()
	IMPLICIT NONE
    
	INTEGER :: ix, iy, iz
    INTEGER(Long) :: thiscell
    REAL(Double) :: aux
    
    U_s = 0
    DO iz = zstart, zend
        DO iy = ystart, yend
            DO ix = xstart, xend
                thiscell=UROLL3(ix, iy, iz)
                S1 = s_E(1,1)*T1(ix, iy, iz)+s_E(1,2)*T2(ix, iy, iz)+s_E(1,3)*T3(ix, iy, iz)
                S2 = s_E(2,1)*T1(ix, iy, iz)+s_E(2,2)*T2(ix, iy, iz)+s_E(2,3)*T3(ix, iy, iz)
                S3 = s_E(3,1)*T1(ix, iy, iz)+s_E(3,2)*T2(ix, iy, iz)+s_E(3,3)*T3(ix, iy, iz)
                S4 = s_E(4,4)*T4(ix, iy, iz)
                S5 = s_E(5,5)*T5(ix, iy, iz)
                S6 = s_E(6,6)*T6(ix, iy, iz)
                aux = S1*(S1*c_E(1,1)+S2*c_E(1,2)+S3*c_E(1,3)) +&
                      S2*(S1*c_E(2,1)+S2*c_E(2,2)+S3*c_E(2,3)) +&
                      S3*(S1*c_E(3,1)+S2*c_E(3,2)+S3*c_E(3,3)) +&
                      S4**2*c_E(4,4) +&
                      S5**2*c_E(5,5) +&
                      S6**2*c_E(6,6)
                U_s=U_s+aux
            END DO
        END DO
    END DO
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(U_s, U_s_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    IF (me .EQ. 0) U_s_total=0.5*U_s_total
    
ENDSUBROUTINE Get_Total_Strain_Energy

SUBROUTINE Get_Total_Electric_Energy()
	IMPLICIT NONE
    
	INTEGER :: ix, iy, iz
    INTEGER(Long) :: thiscell
    
    U_e = 0
    DO iz = zstart, zend
        DO iy = ystart, yend
            DO ix = xstart, xend
                thiscell=UROLL3(ix, iy, iz)
                U_e = U_e + Ex(ix, iy, iz)**2 + Ey(ix, iy, iz)**2 + Ez(ix, iy, iz)**2
            END DO
        END DO
    END DO
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(U_e, U_e_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    IF (me .EQ. 0) U_e_total = 0.5*U_e_total/beta_s(1,1)*dt**2
    
ENDSUBROUTINE Get_Total_Electric_Energy

SUBROUTINE freesurface_to_vtk(vtkfile)
	IMPLICIT NONE
    CHARACTER (LEN=*),INTENT(IN) :: vtkfile
	INTEGER :: ix, iy, iz
    INTEGER :: E_IO
    REAL, DIMENSION(1:Nx-2) :: x_xml_rect 
    REAL, DIMENSION(1:Ny-2) :: y_xml_rect 
    REAL, DIMENSION(1:1) :: z_xml_rect 
    
    ! esempio di output in "RectilinearGrid" XML (binario)
    ! creazione del file
    E_IO = VTK_INI_XML(output_format = 'BINARY',              &
                       filename      = trim(output_dir)//'freesurface'//vtkfile, &
                       mesh_topology = 'RectilinearGrid',     &
                       nx1=1,nx2=nx-2,ny1=1,ny2=ny-2,nz1=1,nz2=1)
    ! salvataggio della geometria del "Piece" corrente
    x_xml_rect=((/(ix,ix=1, Nx-2)/)-1)*deltax+offsetx
    y_xml_rect=((/(iy,iy=1, Ny-2)/)-1)*deltay+offsety
    z_xml_rect=((/(iz,iz=1, 1)/)-1)*deltaz
    E_IO = VTK_GEO_XML(nx1=1,nx2=nx-2,ny1=1,ny2=ny-2,nz1=1,nz2=1, &
                         X=x_xml_rect,Y=y_xml_rect,Z=z_xml_rect)
    ! inizializzazione del salvataggio delle variabili definite ai nodi del "Piece" corrente
    E_IO = VTK_DAT_XML(var_location     = 'node', &
                       var_block_action = 'OPEN')
    ! salvataggio delle variabili definite ai nodi del "Piece" corrente
    E_IO = VTK_VAR_XML(NC_NN   = (Nx-2)*(Ny-2), &
                       varname = 'V',                    &
                       varX=Vx(1:Nx-2, 1:Ny-2, 1:1),varY=Vy(1:Nx-2, 1:Ny-2, 1:1),varZ=Vz(1:Nx-2, 1:Ny-2, 1:1))
    E_IO = VTK_VAR_XML(NC_NN   = (Nx-2)*(Ny-2), &
                       varname = 'T_axial',                    &
                       varX=T1(1:Nx-2, 1:Ny-2, 1:1),varY=T2(1:Nx-2, 1:Ny-2, 1:1),varZ=T3(1:Nx-2, 1:Ny-2, 1:1))
    E_IO = VTK_VAR_XML(NC_NN   = (Nx-2)*(Ny-2), &
                       varname = 'T_cortante',                    &
                       varX=T4(1:Nx-2, 1:Ny-2, 1:1),varY=T5(1:Nx-2, 1:Ny-2, 1:1),varZ=T6(1:Nx-2, 1:Ny-2, 1:1))
    E_IO = VTK_VAR_XML(NC_NN   = (Nx-2)*(Ny-2), &
                       varname = 'E',                    &
                       varX=Ex(1:Nx-2, 1:Ny-2, 1:1),varY=Ey(1:Nx-2, 1:Ny-2, 1:1),varZ=Ez(1:Nx-2, 1:Ny-2, 1:1))
    E_IO = VTK_VAR_XML(NC_NN   = (Nx-2)*(Ny-2), &
                       varname = 'D',                    &
                       varX=Disx(1:Nx-2, 1:Ny-2, 1:1),varY=Disy(1:Nx-2, 1:Ny-2, 1:1),varZ=Disz(1:Nx-2, 1:Ny-2, 1:1))
    ! chiusura del blocco delle variabili definite ai nodi del "Piece" corrente
    E_IO = VTK_DAT_XML(var_location     = 'node', &
                       var_block_action = 'Close')
    ! chiusura del "Piece" corrente
    E_IO = VTK_GEO_XML()
    ! chiusura del file
    E_IO = VTK_END_XML()
    
ENDSUBROUTINE freesurface_to_vtk

SUBROUTINE save_D0_to_vtk(vtkfile, nx1, nx2, ny1, ny2, nz1, nz2)
	IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN) :: vtkfile
	INTEGER, INTENT(IN) :: nx1, nx2, ny1, ny2, nz1, nz2
    
	INTEGER :: ix, iy, iz
    INTEGER :: E_IO
    REAL, DIMENSION(nx1:nx2) :: x_xml_rect 
    REAL, DIMENSION(ny1:ny2) :: y_xml_rect 
    REAL, DIMENSION(nz1:nz2) :: z_xml_rect 
    
    ! esempio di output in "RectilinearGrid" XML (binario)
    ! creazione del file
    E_IO = VTK_INI_XML(output_format = 'BINARY',              &
                       filename      = trim(output_dir)//vtkfile, &
                       mesh_topology = 'RectilinearGrid',     &
                       nx1=nx1,nx2=nx2,ny1=ny1,ny2=ny2,nz1=nz1,nz2=nz2)
    ! salvataggio della geometria del "Piece" corrente
    x_xml_rect=((/(ix,ix=nx1, nx2)/)-1)*deltax+offsetx
    y_xml_rect=((/(iy,iy=ny1, ny2)/)-1)*deltay+offsety
    z_xml_rect=((/(iz,iz=nz1, nz2)/)-1)*deltaz
    E_IO = VTK_GEO_XML(nx1=nx1,nx2=nx2,ny1=ny1,ny2=ny2,nz1=nz1,nz2=nz2, &
                         X=x_xml_rect,Y=y_xml_rect,Z=z_xml_rect)
    ! inizializzazione del salvataggio delle variabili definite ai nodi del "Piece" corrente
    E_IO = VTK_DAT_XML(var_location     = 'node', &
                       var_block_action = 'OPEN')
    ! salvataggio delle variabili definite ai nodi del "Piece" corrente
    E_IO = VTK_VAR_XML(NC_NN   = (nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1), &
                       varname = 'D0',                    &
                       varX=D0x(nx1:nx2,ny1:ny2,nz1:nz2),varY=D0y(nx1:nx2, ny1:ny2, nz1:nz2),varZ=D0z(nx1:nx2, ny1:ny2, nz1:nz2))
    ! chiusura del blocco delle variabili definite ai nodi del "Piece" corrente
    E_IO = VTK_DAT_XML(var_location     = 'node', &
                       var_block_action = 'Close')
    ! chiusura del "Piece" corrente
    E_IO = VTK_GEO_XML()
    ! chiusura del file 
    E_IO = VTK_END_XML()
    
ENDSUBROUTINE save_D0_to_vtk

SUBROUTINE EA_to_vtk(vtkfile, nx1, nx2, ny1, ny2, nz1, nz2)

	IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN) :: vtkfile
	INTEGER, INTENT(IN) :: nx1, nx2, ny1, ny2, nz1, nz2
    
	INTEGER :: ix, iy, iz
    INTEGER :: E_IO
    REAL, DIMENSION(nx1:nx2) :: x_xml_rect 
    REAL, DIMENSION(ny1:ny2) :: y_xml_rect 
    REAL, DIMENSION(nz1:nz2) :: z_xml_rect 
    
    ! esempio di output in "RectilinearGrid" XML (binario)
    ! creazione del file
    E_IO = VTK_INI_XML(output_format = 'BINARY',              &
                       filename      = trim(output_dir)//vtkfile, &
                       mesh_topology = 'RectilinearGrid',     &
                       nx1=nx1,nx2=nx2,ny1=ny1,ny2=ny2,nz1=nz1,nz2=nz2)
    ! salvataggio della geometria del "Piece" corrente
    x_xml_rect=((/(ix,ix=nx1, nx2)/)-1)*deltax+offsetx
    y_xml_rect=((/(iy,iy=ny1, ny2)/)-1)*deltay+offsety
    z_xml_rect=((/(iz,iz=nz1, nz2)/)-1)*deltaz
    E_IO = VTK_GEO_XML(nx1=nx1,nx2=nx2,ny1=ny1,ny2=ny2,nz1=nz1,nz2=nz2, &
                         X=x_xml_rect,Y=y_xml_rect,Z=z_xml_rect)
    ! inizializzazione del salvataggio delle variabili definite ai nodi del "Piece" corrente
    E_IO = VTK_DAT_XML(var_location     = 'node', &
                       var_block_action = 'OPEN')
    ! salvataggio delle variabili definite ai nodi del "Piece" corrente
    E_IO = VTK_VAR_XML(NC_NN   = (nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1), &
                       varname = 'V',                    &
                       varX=Vx(nx1:nx2,ny1:ny2,nz1:nz2),varY=Vy(nx1:nx2, ny1:ny2, nz1:nz2),varZ=Vz(nx1:nx2, ny1:ny2, nz1:nz2))
    E_IO = VTK_VAR_XML(NC_NN   = (nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1), &
                       varname = 'T_axial',                    &
                       varX=T1(nx1:nx2,ny1:ny2,nz1:nz2),varY=T2(nx1:nx2, ny1:ny2, nz1:nz2),varZ=T3(nx1:nx2, ny1:ny2, nz1:nz2))
    E_IO = VTK_VAR_XML(NC_NN   = (nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1), &
                       varname = 'T_cortante',                    &
                       varX=T4(nx1:nx2,ny1:ny2,nz1:nz2),varY=T5(nx1:nx2, ny1:ny2, nz1:nz2),varZ=T6(nx1:nx2, ny1:ny2, nz1:nz2))
    E_IO = VTK_VAR_XML(NC_NN   = (nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1), &
                       varname = 'E',                    &
                       varX=Ex(nx1:nx2,ny1:ny2,nz1:nz2),varY=Ey(nx1:nx2, ny1:ny2, nz1:nz2),varZ=Ez(nx1:nx2, ny1:ny2, nz1:nz2))
    E_IO = VTK_VAR_XML(NC_NN   = (nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1), &
                       varname = 'D',                    &
                       varX=Disx(nx1:nx2,ny1:ny2,nz1:nz2),varY=Disy(nx1:nx2, ny1:ny2, nz1:nz2),varZ=Disz(nx1:nx2, ny1:ny2, nz1:nz2))
    ! chiusura del blocco delle variabili definite ai nodi del "Piece" corrente
    E_IO = VTK_DAT_XML(var_location     = 'node', &
                       var_block_action = 'Close')
    ! chiusura del "Piece" corrente
    E_IO = VTK_GEO_XML()
    ! chiusura del file 
    E_IO = VTK_END_XML()
    
ENDSUBROUTINE EA_to_vtk

SUBROUTINE read_idt(idt_file,ND0x,ND0y,ND0z)
	IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN) :: idt_file
    INTEGER, INTENT(IN) :: ND0x, ND0y, ND0z
    
    INTEGER :: st
    INTEGER :: ix, iy, iz
    INTEGER :: procix, prociy
    REAL(Double), DIMENSION(:,:,:,:), POINTER :: read_D0
    REAL(Double), DIMENSION(3) :: junk
    CHARACTER (LEN=name_len) :: outfile
    
    IF (me==0) THEN
        
    OPEN(UNIT=10,FILE=idt_file,STATUS='old',ACTION='read',IOSTAT=st)
    ALLOCATE(read_D0(0:2, 0:ND0x-1 , 0:ND0y-1 , 0:ND0z-1))
    
    DO iz = 0 , ND0z-1
    DO iy = 0 , ND0y-1
    DO ix = 0 , ND0x-1
        READ(10,*) junk, read_D0(0:2,ix,iy,iz)
    END Do
    END DO
    write (*,*) iz
    END DO
    
    CALL SYSTEM('mkdir -p '//trim(output_dir)//'IDT')
        
    DO prociy=0, Nprocsy-1
    DO procix=0, Nprocsx-1
    
        DO iz=1, Nz-2
        DO iy=1, Ny-2
        DO ix=1, Nx-2
            D0x(ix, iy, iz)=read_D0(0,ix-1+procix*(Nx-2),iy-1+prociy*(Ny-2),iz-1)
            D0y(ix, iy, iz)=read_D0(1,ix-1+procix*(Nx-2),iy-1+prociy*(Ny-2),iz-1)
            D0z(ix, iy, iz)=read_D0(2,ix-1+procix*(Nx-2),iy-1+prociy*(Ny-2),iz-1)
        END DO
        END DO
        END DO
        
        WRITE(outfile, '(A,i3.3)') 'D0',procix+Nprocsx*prociy
        OPEN(UNIT=12, FILE=trim(output_dir)//'IDT/'//outfile, ACTION="write", STATUS="replace", FORM = 'unformatted', IOSTAT = st)
        DO ix = 0, Nx - 1
        WRITE(12) ((D0x(ix, iy, iz), iy = 0, Ny - 1), iz = 0, Nz - 1)
        WRITE(12) ((D0y(ix, iy, iz), iy = 0, Ny - 1), iz = 0, Nz - 1)
        WRITE(12) ((D0z(ix, iy, iz), iy = 0, Ny - 1), iz = 0, Nz - 1)
        END DO
        CLOSE(12)
        
    END DO
    END DO
    
    DEALLOCATE(read_D0,STAT=st)
    
    ENDIF
    
ENDSUBROUTINE read_idt


ENDMODULE Lib_FDTD_SAW
