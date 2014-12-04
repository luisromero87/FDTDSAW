MODULE Lib_FDTD_SAW

USE mpi
USE Type_Kinds
USE Constants_Module
USE Global_Vars

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
PUBLIC :: share_v
PUBLIC :: share_T
PUBLIC :: v_half_step
PUBLIC :: T_half_step
!PUBLIC :: open_vtk_file
!PUBLIC :: write_free_surface
!PUBLIC :: write_volume_v
!PUBLIC :: write_volume_w1
!PUBLIC :: write_volume_D0
!PUBLIC :: xzplane
PUBLIC :: Get_Total_Kinetic_Energy
PUBLIC :: Get_Total_Strain_Energy
PUBLIC :: Get_Total_Electric_Energy

INTERFACE
        SUBROUTINE share_v()
            IMPLICIT NONE
        END SUBROUTINE share_v
        
        SUBROUTINE share_T()
            IMPLICIT NONE
        END SUBROUTINE share_T
        
        SUBROUTINE v_half_step(zper)
!            USE Type_Kinds
!            USE Constants_Module
!            USE Global_Vars
            IMPLICIT NONE
            CHARACTER(LEN=*), OPTIONAL :: zper
        ENDSUBROUTINE v_half_step
        
        SUBROUTINE T_half_step(zper)
!            USE Type_Kinds
!            USE Constants_Module
!            USE Global_Vars
            IMPLICIT NONE
            CHARACTER(LEN=*), OPTIONAL :: zper
        ENDSUBROUTINE T_half_step
ENDINTERFACE

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

SUBROUTINE read_input_param(input_param)
    IMPLICIT NONE

    CHARACTER(LEN = name_len), OPTIONAL :: input_param
    
    INTEGER :: st
    CHARACTER (LEN = name_len) :: param = ''
    
    IF (PRESENT(input_param)) THEN
        sim_config=input_param
    END IF
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
    DO WHILE (param /= 'Data_fstep:')
        READ(12, *) param, data_fstep
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
    
    DO WHILE (param /= 'PML_width:')
        READ(12, *) param, PMLwidth
    END DO
    DO WHILE (param /= 'smax:')
        READ(12, *) param, smax
    END DO
    DO WHILE (param /= 'm:')
        READ(12, *) param, m
    END DO

    CLOSE(12)
    WRITE(*,*) "\n\n"
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
            CALL SYSTEM('mkdir -p outputdata/weights')
            WRITE(which, '(A,I3.3)') 'outputdata/weights/weights', me
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
    
    WRITE(which, '(A,I3.3)') 'outputdata/IDT/D0', me
    OPEN(UNIT = 12, FILE = which, ACTION="read", STATUS="old", FORM = 'unformatted', IOSTAT = st)
    
    IF (me .EQ. 0 .AND. st .EQ. 0) THEN
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

ENDMODULE Lib_FDTD_SAW
