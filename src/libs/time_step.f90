!     
! File:   time_step.f90
! Author: ludwig
!
! Created on April 12, 2014, 9:35 PM
!

SUBROUTINE v_half_step()
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE

    INTEGER(Long) :: UROLL3

    INTEGER :: ix, iy, iz


    DO iz = 1, Nz - 2
        DO iy = 1, Ny - 2
            DO ix = 1, Nx - 2
                Vx_x(UROLL3(ix, iy, iz)) = Vx_x(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), x) + &
                rho_inv * w2(UROLL3(ix, iy, iz), x) * (T1(UROLL3(ix + 1, iy, iz)) - T1(UROLL3(ix, iy, iz)))

                Vx_y(UROLL3(ix, iy, iz)) = Vx_y(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), y) + &
                rho_inv * w2(UROLL3(ix, iy, iz), y) * (T6(UROLL3(ix, iy, iz)) - T6(UROLL3(ix, iy - 1, iz)))

                Vx_z(UROLL3(ix, iy, iz)) = Vx_z(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), z) + &
                rho_inv * w2(UROLL3(ix, iy, iz), z) * (T5(UROLL3(ix, iy, iz)) - T5(UROLL3(ix, iy, iz - 1)))

                Vx(UROLL3(ix, iy, iz)) = Vx_x(UROLL3(ix, iy, iz)) + Vx_y(UROLL3(ix, iy, iz)) + Vx_z(UROLL3(ix, iy, iz))


                Vy_x(UROLL3(ix, iy, iz)) = Vy_x(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), x) + &
                rho_inv * w2(UROLL3(ix, iy, iz), x) * (T6(UROLL3(ix, iy, iz)) - T6(UROLL3(ix - 1, iy, iz)))

                Vy_y(UROLL3(ix, iy, iz)) = Vy_y(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), y) + &
                rho_inv * w2(UROLL3(ix, iy, iz), y) * (T2(UROLL3(ix, iy + 1, iz)) - T2(UROLL3(ix, iy, iz)))

                Vy_z(UROLL3(ix, iy, iz)) = Vy_z(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), z) + &
                rho_inv * w2(UROLL3(ix, iy, iz), z) * (T4(UROLL3(ix, iy, iz)) - T4(UROLL3(ix, iy, iz - 1)))

                Vy(UROLL3(ix, iy, iz)) = Vy_x(UROLL3(ix, iy, iz)) + Vy_y(UROLL3(ix, iy, iz)) + Vy_z(UROLL3(ix, iy, iz))


                Vz_x(UROLL3(ix, iy, iz)) = Vz_x(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), x) + &
                rho_inv * w2(UROLL3(ix, iy, iz), x) * (T5(UROLL3(ix, iy, iz)) - T5(UROLL3(ix - 1, iy, iz)))

                Vz_y(UROLL3(ix, iy, iz)) = Vz_y(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), y) + &
                rho_inv * w2(UROLL3(ix, iy, iz), y) * (T4(UROLL3(ix, iy, iz)) - T4(UROLL3(ix, iy - 1, iz)))

                Vz_z(UROLL3(ix, iy, iz)) = Vz_z(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), z) + &
                rho_inv * w2(UROLL3(ix, iy, iz), z) * (T3(UROLL3(ix, iy, iz + 1)) - T3(UROLL3(ix, iy, iz)))

                Vz(UROLL3(ix, iy, iz)) = Vz_x(UROLL3(ix, iy, iz)) + Vz_y(UROLL3(ix, iy, iz)) + Vz_z(UROLL3(ix, iy, iz))
            END DO
        END DO
    END DO


END SUBROUTINE v_half_step


SUBROUTINE free_boundary_v()
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE

    INTEGER(Long) :: UROLL3

    INTEGER :: ix, iy, iz


    iz = 1
    DO iy = 1, Ny - 2
        DO ix = 1, Nx - 2
            Vz(UROLL3(ix, iy, iz - 1)) = Vz(UROLL3(ix, iy, iz))+(dz(iz)/c_E(12 + 3))*(&
            c_E(12 + 1)/dx(ix) * (Vx(UROLL3(ix, iy, iz)) - Vx(UROLL3(ix - 1, iy, iz)))+&
            c_E(12 + 2)/dy(iy) * (Vy(UROLL3(ix, iy, iz)) - Vy(UROLL3(ix, iy - 1, iz))))
        END DO
    END DO
    
END SUBROUTINE free_boundary_v


SUBROUTINE T_half_step()
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE

    INTEGER(Long) :: UROLL3

    INTEGER :: ix, iy, iz

    REAL(Double) :: dVxdx = 0.0_dp
    REAL(Double) :: dVydy = 0.0_dp
    REAL(Double) :: dVzdz = 0.0_dp



    DO iz = 1, Nz - 2
        DO iy = 1, Ny - 2
            DO ix = 1, Nx - 2

                dVxdx = Vx(UROLL3(ix, iy, iz)) - Vx(UROLL3(ix - 1, iy, iz))
                dVydy = Vy(UROLL3(ix, iy, iz)) - Vy(UROLL3(ix, iy - 1, iz))
                dVzdz = Vz(UROLL3(ix, iy, iz)) - Vz(UROLL3(ix, iy, iz - 1))

                T1_x(UROLL3(ix, iy, iz)) = T1_x(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), x) +&
                w2(UROLL3(ix, iy, iz), x) * c_E(1)*(dVxdx)
                T1_y(UROLL3(ix, iy, iz)) = T1_y(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), y) +&
                w2(UROLL3(ix, iy, iz), y) * c_E(2)*(dVydy)
                T1_z(UROLL3(ix, iy, iz)) = T1_z(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), z) +&
                w2(UROLL3(ix, iy, iz), z) * c_E(3)*(dVzdz)
                T1(UROLL3(ix, iy, iz)) = T1_x(UROLL3(ix, iy, iz)) + T1_y(UROLL3(ix, iy, iz)) + T1_z(UROLL3(ix, iy, iz))


                T2_x(UROLL3(ix, iy, iz)) = T2_x(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), x) +&
                w2(UROLL3(ix, iy, iz), x) * c_E(6 + 1)*(dVxdx)
                T2_y(UROLL3(ix, iy, iz)) = T2_y(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), y) +&
                w2(UROLL3(ix, iy, iz), y) * c_E(6 + 2)*(dVydy)
                T2_z(UROLL3(ix, iy, iz)) = T2_z(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), z) +&
                w2(UROLL3(ix, iy, iz), z) * c_E(6 + 3)*(dVzdz)
                T2(UROLL3(ix, iy, iz)) = T2_x(UROLL3(ix, iy, iz)) + T2_y(UROLL3(ix, iy, iz)) + T2_z(UROLL3(ix, iy, iz))

                T3_x(UROLL3(ix, iy, iz)) = T3_x(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), x) +&
                w2(UROLL3(ix, iy, iz), x) * c_E(12 + 1)*(dVxdx)
                T3_y(UROLL3(ix, iy, iz)) = T3_y(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), y) +&
                w2(UROLL3(ix, iy, iz), y) * c_E(12 + 2)*(dVydy)
                T3_z(UROLL3(ix, iy, iz)) = T3_z(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), z) +&
                w2(UROLL3(ix, iy, iz), z) * c_E(12 + 3)*(dVzdz)
                T3(UROLL3(ix, iy, iz)) = T3_x(UROLL3(ix, iy, iz)) + T3_y(UROLL3(ix, iy, iz)) + T3_z(UROLL3(ix, iy, iz))

                T4_y(UROLL3(ix, iy, iz)) = T4_y(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), y) +&
                w2(UROLL3(ix, iy, iz), y) * c_E(18 + 4)*(Vz(UROLL3(ix, iy + 1, iz)) - Vz(UROLL3(ix, iy, iz)))
                T4_z(UROLL3(ix, iy, iz)) = T4_z(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), z) +&
                w2(UROLL3(ix, iy, iz), z) * c_E(18 + 4)*(Vy(UROLL3(ix, iy, iz + 1)) - Vy(UROLL3(ix, iy, iz)))
                T4(UROLL3(ix, iy, iz)) = T4_y(UROLL3(ix, iy, iz)) + T4_z(UROLL3(ix, iy, iz))

                T5_x(UROLL3(ix, iy, iz)) = T5_x(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), x) +&
                w2(UROLL3(ix, iy, iz), x) * c_E(24 + 5)*(Vz(UROLL3(ix + 1, iy, iz)) - Vz(UROLL3(ix, iy, iz)))
                T5_z(UROLL3(ix, iy, iz)) = T5_z(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), z) +&
                w2(UROLL3(ix, iy, iz), z) * c_E(24 + 5)*(Vx(UROLL3(ix, iy, iz + 1)) - Vx(UROLL3(ix, iy, iz)))
                T5(UROLL3(ix, iy, iz)) = T5_x(UROLL3(ix, iy, iz)) + T5_z(UROLL3(ix, iy, iz))

                T6_x(UROLL3(ix, iy, iz)) = T6_x(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), x) +&
                w2(UROLL3(ix, iy, iz), x) * c_E(30 + 6)*(Vy(UROLL3(ix + 1, iy, iz)) - Vy(UROLL3(ix, iy, iz)))
                T6_y(UROLL3(ix, iy, iz)) = T6_y(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), y) +&
                w2(UROLL3(ix, iy, iz), y) * c_E(30 + 6)*(Vx(UROLL3(ix, iy + 1, iz)) - Vx(UROLL3(ix, iy, iz)))
                T6(UROLL3(ix, iy, iz)) = T6_x(UROLL3(ix, iy, iz)) + T6_y(UROLL3(ix, iy, iz))

            END DO
        END DO
    END DO



END SUBROUTINE T_half_step


SUBROUTINE free_boundary_T()
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE

    INTEGER(Long) :: UROLL3

    INTEGER :: ix, iy, iz


    iz = 1
    DO iy = 1, Ny - 2
        DO ix = 1, Nx - 2
            T5(UROLL3(ix, iy, iz - 1)) = -T5(UROLL3(ix, iy, iz))
            T3(UROLL3(ix, iy, iz)) = 0
            T4(UROLL3(ix, iy, iz - 1)) = -T4(UROLL3(ix, iy, iz))
        END DO
    END DO


END SUBROUTINE free_boundary_T

SUBROUTINE dot_source()
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE

    INTEGER(Long) :: UROLL3

    INTEGER :: ix, iy, iz

    ix = Nx-2
    iy = Ny-2
    iz = 1

    Vz(UROLL3(ix, iy, iz)) = Vz(UROLL3(ix, iy, iz)) + &
    (step * dt - 3 * PWIDTH)/(3 * PWIDTH) * exp(-1.0 * ((step * dt - 3.0 * PWIDTH)/(PWIDTH))**2)

    !WRITE(*, *) (step * dt - 3 * PWIDTH)/(3 * PWIDTH) * exp(-1.0 * ((step * dt - 3.0 * PWIDTH)/(PWIDTH))**2)

END SUBROUTINE dot_source




SUBROUTINE share_v()
	USE mpi
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
	IMPLICIT NONE
	INTEGER(Long) :: UROLL3
	INTEGER(Short) UROLLPROC
	INTEGER :: ix, iy, iz, nextprocid
	REAL(Double), DIMENSION (0:3*Nz*Ny-1) :: mpibufferx
	REAL(Double), DIMENSION (0:3*Nz*Nx-1) :: mpibuffery
	
	IF (ntasks>1 .AND. me<Nprocsx*Nprocsy) THEN
		!direccion x
		!interfaces pares
        IF ( MOD(procsx,2)==0 .AND. (procsx<Nprocsx-1) ) THEN
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(0*Nz*Ny+iz*Ny+iy)= Vx(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(1*Nz*Ny+iz*Ny+iy)= Vy(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(2*Nz*Ny+iz*Ny+iy)= Vz(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			CALL MPI_SEND(mpibufferx,3*Nz*Ny,MPI_DOUBLE_PRECISION,int(UROLLPROC(procsx+1,procsy)),me,MPI_COMM_WORLD,ierr)
			
			CALL MPI_RECV(mpibufferx,3*Nz*Ny,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				Vx(UROLL3(Nx-1,iy,iz))=mpibufferx(0*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				Vy(UROLL3(Nx-1,iy,iz))=mpibufferx(1*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				Vz(UROLL3(Nx-1,iy,iz))=mpibufferx(2*Nz*Ny+iz*Ny+iy)
			END DO
			END DO

		ELSE IF ( MOD(procsx,2)==1 ) THEN
		
			CALL MPI_RECV(mpibufferx,3*Nz*Ny,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				Vx(UROLL3(0_li,iy,iz))=mpibufferx(0*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				Vy(UROLL3(0_li,iy,iz))=mpibufferx(1*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				Vz(UROLL3(0_li,iy,iz))=mpibufferx(2*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(0*Nz*Ny+iz*Ny+iy)= Vx(UROLL3(1_li,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(1*Nz*Ny+iz*Ny+iy)= Vy(UROLL3(1_li,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(2*Nz*Ny+iz*Ny+iy)= Vz(UROLL3(1_li,iy,iz))
			END DO
			END DO
			CALL MPI_SEND(mpibufferx,3*Nz*Ny,MPI_DOUBLE_PRECISION,int(UROLLPROC(procsx-1,procsy)),me,MPI_COMM_WORLD,ierr)
			
		END IF
		!interfaces impares
        IF ( MOD(procsx,2)==1 .AND. (procsx<Nprocsx-1) ) THEN
			
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(0*Nz*Ny+iz*Ny+iy)= Vx(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(1*Nz*Ny+iz*Ny+iy)= Vy(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(2*Nz*Ny+iz*Ny+iy)= Vz(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			CALL MPI_SEND(mpibufferx,3*Nz*Ny,MPI_DOUBLE_PRECISION,int(UROLLPROC(procsx+1,procsy)),me,MPI_COMM_WORLD,ierr)
			
			CALL MPI_RECV(mpibufferx,3*Nz*Ny,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				Vx(UROLL3(Nx-1,iy,iz))=mpibufferx(0*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				Vy(UROLL3(Nx-1,iy,iz))=mpibufferx(1*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				Vz(UROLL3(Nx-1,iy,iz))=mpibufferx(2*Nz*Ny+iz*Ny+iy)
			END DO
			END DO

		ELSE IF ( MOD(procsx,2)==0 .AND. (procsx >0) ) THEN
		
			CALL MPI_RECV(mpibufferx,3*Nz*Ny,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				Vx(UROLL3(0_li,iy,iz))=mpibufferx(0*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				Vy(UROLL3(0_li,iy,iz))=mpibufferx(1*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				Vz(UROLL3(0_li,iy,iz))=mpibufferx(2*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(0*Nz*Ny+iz*Ny+iy)= Vx(UROLL3(1_li,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(1*Nz*Ny+iz*Ny+iy)= Vy(UROLL3(1_li,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(2*Nz*Ny+iz*Ny+iy)= Vz(UROLL3(1_li,iy,iz))
			END DO
			END DO
			CALL MPI_SEND(mpibufferx,3*Nz*Ny,MPI_DOUBLE_PRECISION,int(UROLLPROC(procsx-1,procsy)),me,MPI_COMM_WORLD,ierr)
			
		END IF
		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
		!direccion y
		!interfaces pares
        IF ( MOD(procsy,2)==0 .AND. (procsy<Nprocsy-1) ) THEN
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(0*Nz*Nx+iz*Nx+ix)= Vx(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(1*Nz*Nx+iz*Nx+ix)= Vy(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(2*Nz*Nx+iz*Nx+ix)= Vz(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			CALL MPI_SEND(mpibuffery,3*Nz*Nx,MPI_DOUBLE_PRECISION,int(UROLLPROC(procsx,procsy+1)),me,MPI_COMM_WORLD,ierr)
			
			CALL MPI_RECV(mpibuffery,3*Nz*Nx,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				Vx(UROLL3(ix,Ny-1,iz))=mpibuffery(0*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				Vy(UROLL3(ix,Ny-1,iz))=mpibuffery(1*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				Vz(UROLL3(ix,Ny-1,iz))=mpibuffery(2*Nz*Nx+iz*Nx+ix)
			END DO
			END DO

		ELSE IF ( MOD(procsy,2)==1 ) THEN
		
			CALL MPI_RECV(mpibuffery,3*Nz*Nx,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				Vx(UROLL3(ix,0_li,iz))=mpibuffery(0*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				Vy(UROLL3(ix,0_li,iz))=mpibuffery(1*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				Vz(UROLL3(ix,0_li,iz))=mpibuffery(2*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(0*Nz*Nx+iz*Nx+ix)= Vx(UROLL3(ix,1_li,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(1*Nz*Nx+iz*Nx+ix)= Vy(UROLL3(ix,1_li,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(2*Nz*Nx+iz*Nx+ix)= Vz(UROLL3(ix,1_li,iz))
			END DO
			END DO
			CALL MPI_SEND(mpibuffery,3*Nz*Nx,MPI_DOUBLE_PRECISION,int(UROLLPROC(procsx,procsy-1)),me,MPI_COMM_WORLD,ierr)
			
		END IF
		!interfaces impares
        IF ( MOD(procsy,2)==1 .AND. (procsy<Nprocsy-1) ) THEN
			
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(0*Nz*Nx+iz*Nx+ix)= Vx(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(1*Nz*Nx+iz*Nx+ix)= Vy(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(2*Nz*Nx+iz*Nx+ix)= Vz(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			CALL MPI_SEND(mpibuffery,3*Nz*Nx,MPI_DOUBLE_PRECISION,int(UROLLPROC(procsx,procsy+1)),me,MPI_COMM_WORLD,ierr)
			
			CALL MPI_RECV(mpibuffery,3*Nz*Nx,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				Vx(UROLL3(ix,Ny-1,iz))=mpibuffery(0*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				Vy(UROLL3(ix,Ny-1,iz))=mpibuffery(1*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				Vz(UROLL3(ix,Ny-1,iz))=mpibuffery(2*Nz*Nx+iz*Nx+ix)
			END DO
			END DO

		ELSE IF ( MOD(procsy,2)==0 .AND. (procsy >0) ) THEN
		
			CALL MPI_RECV(mpibuffery,3*Nz*Nx,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				Vx(UROLL3(ix,0_li,iz))=mpibuffery(0*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				Vy(UROLL3(ix,0_li,iz))=mpibuffery(1*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				Vz(UROLL3(ix,0_li,iz))=mpibuffery(2*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(0*Nz*Nx+iz*Nx+ix)= Vx(UROLL3(ix,1_li,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(1*Nz*Nx+iz*Nx+ix)= Vy(UROLL3(ix,1_li,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(2*Nz*Nx+iz*Nx+ix)= Vz(UROLL3(ix,1_li,iz))
			END DO
			END DO
			CALL MPI_SEND(mpibuffery,3*Nz*Nx,MPI_DOUBLE_PRECISION,int(UROLLPROC(procsx,procsy-1)),me,MPI_COMM_WORLD,ierr)
			
		END IF
	END IF

END SUBROUTINE share_v


SUBROUTINE share_T()
	USE mpi
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
	IMPLICIT NONE
	INTEGER(Long) :: UROLL3
	INTEGER(Short) UROLLPROC
	INTEGER :: ix, iy, iz, nextprocid
	REAL(Double), DIMENSION (0:6*Nz*Ny-1) :: mpibufferx
	REAL(Double), DIMENSION (0:6*Nz*Nx-1) :: mpibuffery
	
	IF (ntasks>1 .AND. me<Nprocsx*Nprocsy) THEN
		!direccion x
		!interfaces pares
        IF ( MOD(procsx,2)==0 .AND. (procsx<Nprocsx-1) ) THEN
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(0*Nz*Ny+iz*Ny+iy)= T1(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(1*Nz*Ny+iz*Ny+iy)= T2(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(2*Nz*Ny+iz*Ny+iy)= T3(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(3*Nz*Ny+iz*Ny+iy)= T4(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(4*Nz*Ny+iz*Ny+iy)= T5(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(5*Nz*Ny+iz*Ny+iy)= T6(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			CALL MPI_SEND(mpibufferx,6*Nz*Ny,MPI_DOUBLE_PRECISION,int(UROLLPROC(procsx+1,procsy)),me,MPI_COMM_WORLD,ierr)
			
			CALL MPI_RECV(mpibufferx,6*Nz*Ny,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T1(UROLL3(Nx-1,iy,iz))=mpibufferx(0*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T2(UROLL3(Nx-1,iy,iz))=mpibufferx(1*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T3(UROLL3(Nx-1,iy,iz))=mpibufferx(2*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T4(UROLL3(Nx-1,iy,iz))=mpibufferx(3*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T5(UROLL3(Nx-1,iy,iz))=mpibufferx(4*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T6(UROLL3(Nx-1,iy,iz))=mpibufferx(5*Nz*Ny+iz*Ny+iy)
			END DO
			END DO

		ELSE IF ( MOD(procsx,2)==1 ) THEN
		
			CALL MPI_RECV(mpibufferx,6*Nz*Ny,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T1(UROLL3(0_li,iy,iz))=mpibufferx(0*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T2(UROLL3(0_li,iy,iz))=mpibufferx(1*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T3(UROLL3(0_li,iy,iz))=mpibufferx(2*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T4(UROLL3(0_li,iy,iz))=mpibufferx(3*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T5(UROLL3(0_li,iy,iz))=mpibufferx(4*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T6(UROLL3(0_li,iy,iz))=mpibufferx(5*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(0*Nz*Ny+iz*Ny+iy)= T1(UROLL3(1_li,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(1*Nz*Ny+iz*Ny+iy)= T2(UROLL3(1_li,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(2*Nz*Ny+iz*Ny+iy)= T3(UROLL3(1_li,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(3*Nz*Ny+iz*Ny+iy)= T4(UROLL3(1_li,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(4*Nz*Ny+iz*Ny+iy)= T5(UROLL3(1_li,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(5*Nz*Ny+iz*Ny+iy)= T6(UROLL3(1_li,iy,iz))
			END DO
			END DO
			CALL MPI_SEND(mpibufferx,6*Nz*Ny,MPI_DOUBLE_PRECISION,int(UROLLPROC(procsx-1,procsy)),me,MPI_COMM_WORLD,ierr)
			
		END IF
		!interfaces impares
        IF ( MOD(procsx,2)==1 .AND. (procsx<Nprocsx-1) ) THEN
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(0*Nz*Ny+iz*Ny+iy)= T1(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(1*Nz*Ny+iz*Ny+iy)= T2(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(2*Nz*Ny+iz*Ny+iy)= T3(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(3*Nz*Ny+iz*Ny+iy)= T4(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(4*Nz*Ny+iz*Ny+iy)= T5(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(5*Nz*Ny+iz*Ny+iy)= T6(UROLL3(Nx-2,iy,iz))
			END DO
			END DO
			CALL MPI_SEND(mpibufferx,6*Nz*Ny,MPI_DOUBLE_PRECISION,int(UROLLPROC(procsx+1,procsy)),me,MPI_COMM_WORLD,ierr)
			
			CALL MPI_RECV(mpibufferx,6*Nz*Ny,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T1(UROLL3(Nx-1,iy,iz))=mpibufferx(0*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T2(UROLL3(Nx-1,iy,iz))=mpibufferx(1*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T3(UROLL3(Nx-1,iy,iz))=mpibufferx(2*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T4(UROLL3(Nx-1,iy,iz))=mpibufferx(3*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T5(UROLL3(Nx-1,iy,iz))=mpibufferx(4*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T6(UROLL3(Nx-1,iy,iz))=mpibufferx(5*Nz*Ny+iz*Ny+iy)
			END DO
			END DO

		ELSE IF ( MOD(procsx,2)==0 .AND. (procsx >0) ) THEN
		
			CALL MPI_RECV(mpibufferx,6*Nz*Ny,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T1(UROLL3(0_li,iy,iz))=mpibufferx(0*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T2(UROLL3(0_li,iy,iz))=mpibufferx(1*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T3(UROLL3(0_li,iy,iz))=mpibufferx(2*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T4(UROLL3(0_li,iy,iz))=mpibufferx(3*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T5(UROLL3(0_li,iy,iz))=mpibufferx(4*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
				T6(UROLL3(0_li,iy,iz))=mpibufferx(5*Nz*Ny+iz*Ny+iy)
			END DO
			END DO
			
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(0*Nz*Ny+iz*Ny+iy)= T1(UROLL3(1_li,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(1*Nz*Ny+iz*Ny+iy)= T2(UROLL3(1_li,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(2*Nz*Ny+iz*Ny+iy)= T3(UROLL3(1_li,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(3*Nz*Ny+iz*Ny+iy)= T4(UROLL3(1_li,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(4*Nz*Ny+iz*Ny+iy)= T5(UROLL3(1_li,iy,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO iy=0, Ny-1
					mpibufferx(5*Nz*Ny+iz*Ny+iy)= T6(UROLL3(1_li,iy,iz))
			END DO
			END DO
			CALL MPI_SEND(mpibufferx,6*Nz*Ny,MPI_DOUBLE_PRECISION,int(UROLLPROC(procsx-1,procsy)),me,MPI_COMM_WORLD,ierr)
			
		END IF
		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
		!direccion y
		!interfaces pares
        IF ( MOD(procsy,2)==0 .AND. (procsy<Nprocsy-1) ) THEN
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(0*Nz*Nx+iz*Nx+ix)= T1(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(1*Nz*Nx+iz*Nx+ix)= T2(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(2*Nz*Nx+iz*Nx+ix)= T3(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(3*Nz*Nx+iz*Nx+ix)= T4(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(4*Nz*Nx+iz*Nx+ix)= T5(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(5*Nz*Nx+iz*Nx+ix)= T6(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			CALL MPI_SEND(mpibuffery,6*Nz*Nx,MPI_DOUBLE_PRECISION,int(UROLLPROC(procsx,procsy+1)),me,MPI_COMM_WORLD,ierr)
			
			CALL MPI_RECV(mpibuffery,6*Nz*Nx,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T1(UROLL3(ix,Ny-1,iz))=mpibuffery(0*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T2(UROLL3(ix,Ny-1,iz))=mpibuffery(1*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T3(UROLL3(ix,Ny-1,iz))=mpibuffery(2*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T4(UROLL3(ix,Ny-1,iz))=mpibuffery(3*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T5(UROLL3(ix,Ny-1,iz))=mpibuffery(4*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T6(UROLL3(ix,Ny-1,iz))=mpibuffery(5*Nz*Nx+iz*Nx+ix)
			END DO
			END DO

		ELSE IF ( MOD(procsy,2)==1 ) THEN
		
			CALL MPI_RECV(mpibuffery,6*Nz*Nx,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T1(UROLL3(ix,0_li,iz))=mpibuffery(0*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T2(UROLL3(ix,0_li,iz))=mpibuffery(1*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T3(UROLL3(ix,0_li,iz))=mpibuffery(2*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T4(UROLL3(ix,0_li,iz))=mpibuffery(3*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T5(UROLL3(ix,0_li,iz))=mpibuffery(4*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T6(UROLL3(ix,0_li,iz))=mpibuffery(5*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(0*Nz*Nx+iz*Nx+ix)= T1(UROLL3(ix,1_li,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(1*Nz*Nx+iz*Nx+ix)= T2(UROLL3(ix,1_li,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(2*Nz*Nx+iz*Nx+ix)= T3(UROLL3(ix,1_li,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(3*Nz*Nx+iz*Nx+ix)= T4(UROLL3(ix,1_li,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(4*Nz*Nx+iz*Nx+ix)= T5(UROLL3(ix,1_li,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(5*Nz*Nx+iz*Nx+ix)= T6(UROLL3(ix,1_li,iz))
			END DO
			END DO
			CALL MPI_SEND(mpibuffery,6*Nz*Nx,MPI_DOUBLE_PRECISION,int(UROLLPROC(procsx,procsy-1)),me,MPI_COMM_WORLD,ierr)
			
		END IF
		!interfaces impares
        IF ( MOD(procsy,2)==1 .AND. (procsy<Nprocsy-1) ) THEN
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(0*Nz*Nx+iz*Nx+ix)= T1(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(1*Nz*Nx+iz*Nx+ix)= T2(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(2*Nz*Nx+iz*Nx+ix)= T3(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(3*Nz*Nx+iz*Nx+ix)= T4(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(4*Nz*Nx+iz*Nx+ix)= T5(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(5*Nz*Nx+iz*Nx+ix)= T6(UROLL3(ix,Ny-2,iz))
			END DO
			END DO
			CALL MPI_SEND(mpibuffery,6*Nz*Nx,MPI_DOUBLE_PRECISION,int(UROLLPROC(procsx,procsy+1)),me,MPI_COMM_WORLD,ierr)
			
			CALL MPI_RECV(mpibuffery,6*Nz*Nx,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T1(UROLL3(ix,Ny-1,iz))=mpibuffery(0*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T2(UROLL3(ix,Ny-1,iz))=mpibuffery(1*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T3(UROLL3(ix,Ny-1,iz))=mpibuffery(2*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T4(UROLL3(ix,Ny-1,iz))=mpibuffery(3*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T5(UROLL3(ix,Ny-1,iz))=mpibuffery(4*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T6(UROLL3(ix,Ny-1,iz))=mpibuffery(5*Nz*Nx+iz*Nx+ix)
			END DO
			END DO

		ELSE IF ( MOD(procsy,2)==0 .AND. (procsy >0) ) THEN
		
			CALL MPI_RECV(mpibuffery,6*Nz*Nx,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,mpistatus,ierr)
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T1(UROLL3(ix,0_li,iz))=mpibuffery(0*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T2(UROLL3(ix,0_li,iz))=mpibuffery(1*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T3(UROLL3(ix,0_li,iz))=mpibuffery(2*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T4(UROLL3(ix,0_li,iz))=mpibuffery(3*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T5(UROLL3(ix,0_li,iz))=mpibuffery(4*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
				T6(UROLL3(ix,0_li,iz))=mpibuffery(5*Nz*Nx+iz*Nx+ix)
			END DO
			END DO
			
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(0*Nz*Nx+iz*Nx+ix)= T1(UROLL3(ix,1_li,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(1*Nz*Nx+iz*Nx+ix)= T2(UROLL3(ix,1_li,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(2*Nz*Nx+iz*Nx+ix)= T3(UROLL3(ix,1_li,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(3*Nz*Nx+iz*Nx+ix)= T4(UROLL3(ix,1_li,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(4*Nz*Nx+iz*Nx+ix)= T5(UROLL3(ix,1_li,iz))
			END DO
			END DO
			DO iz=0, Nz-1
			DO ix=0, Nx-1
					mpibuffery(5*Nz*Nx+iz*Nx+ix)= T6(UROLL3(ix,1_li,iz))
			END DO
			END DO
			CALL MPI_SEND(mpibuffery,6*Nz*Nx,MPI_DOUBLE_PRECISION,int(UROLLPROC(procsx,procsy-1)),me,MPI_COMM_WORLD,ierr)
			
		END IF
	END IF

END SUBROUTINE share_T









