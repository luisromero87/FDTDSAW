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
    INTEGER(Long) :: thiscell, anothercell
    REAL (Double) :: w1total
    
    
    DO iy=0, Ny-1
	DO ix=0, Nx-1
		iz=0
		thiscell=UROLL3(ix, iy, iz)
		iz=Nz-2
		anothercell=UROLL3(ix, iy, iz)
		T1(thiscell)=T1(anothercell)
		T2(thiscell)=T2(anothercell)
		T3(thiscell)=T3(anothercell)
		T4(thiscell)=T4(anothercell)
		T5(thiscell)=T5(anothercell)
		T6(thiscell)=T6(anothercell)
		iz=Nz-1
		thiscell=UROLL3(ix, iy, iz)
		iz=1
		anothercell=UROLL3(ix, iy, iz)
		T1(thiscell)=T1(anothercell)
		T2(thiscell)=T2(anothercell)
		T3(thiscell)=T3(anothercell)
		T4(thiscell)=T4(anothercell)
		T5(thiscell)=T5(anothercell)
		T6(thiscell)=T6(anothercell)
    END DO
    END DO
        
    CALL MPI_WAITALL(4, reqs, stats, ierr)
    
    DO iz=0, Nz-1
    DO iy=0, Ny-1
        T1(UROLL3(Nx-1,iy,iz))=recvbuff_RIGHT(0*Nz*Ny+iz*Ny+iy)
        T2(UROLL3(Nx-1,iy,iz))=recvbuff_RIGHT(1*Nz*Ny+iz*Ny+iy)
        T3(UROLL3(Nx-1,iy,iz))=recvbuff_RIGHT(2*Nz*Ny+iz*Ny+iy)
        T4(UROLL3(Nx-1,iy,iz))=recvbuff_RIGHT(3*Nz*Ny+iz*Ny+iy)
        T5(UROLL3(Nx-1,iy,iz))=recvbuff_RIGHT(4*Nz*Ny+iz*Ny+iy)
        T6(UROLL3(Nx-1,iy,iz))=recvbuff_RIGHT(5*Nz*Ny+iz*Ny+iy)
        T1(UROLL3(0_li,iy,iz))= recvbuff_LEFT(0*Nz*Ny+iz*Ny+iy)
        T2(UROLL3(0_li,iy,iz))= recvbuff_LEFT(1*Nz*Ny+iz*Ny+iy)
        T3(UROLL3(0_li,iy,iz))= recvbuff_LEFT(2*Nz*Ny+iz*Ny+iy)
        T4(UROLL3(0_li,iy,iz))= recvbuff_LEFT(3*Nz*Ny+iz*Ny+iy)
        T5(UROLL3(0_li,iy,iz))= recvbuff_LEFT(4*Nz*Ny+iz*Ny+iy)
        T6(UROLL3(0_li,iy,iz))= recvbuff_LEFT(5*Nz*Ny+iz*Ny+iy)
    END DO
    END DO
    
    DO iz=0, Nz-1
    DO ix=0, Nx-1
        T1(UROLL3(ix,Ny-1,iz))=recvbuff_UP(0*Nz*Nx+iz*Nx+ix)
        T2(UROLL3(ix,Ny-1,iz))=recvbuff_UP(1*Nz*Nx+iz*Nx+ix)
        T3(UROLL3(ix,Ny-1,iz))=recvbuff_UP(2*Nz*Nx+iz*Nx+ix)
        T4(UROLL3(ix,Ny-1,iz))=recvbuff_UP(3*Nz*Nx+iz*Nx+ix)
        T5(UROLL3(ix,Ny-1,iz))=recvbuff_UP(4*Nz*Nx+iz*Nx+ix)
        T6(UROLL3(ix,Ny-1,iz))=recvbuff_UP(5*Nz*Nx+iz*Nx+ix)
        T1(UROLL3(ix,0_li,iz))= recvbuff_DOWN(0*Nz*Nx+iz*Nx+ix)
        T2(UROLL3(ix,0_li,iz))= recvbuff_DOWN(1*Nz*Nx+iz*Nx+ix)
        T3(UROLL3(ix,0_li,iz))= recvbuff_DOWN(2*Nz*Nx+iz*Nx+ix)
        T4(UROLL3(ix,0_li,iz))= recvbuff_DOWN(3*Nz*Nx+iz*Nx+ix)
        T5(UROLL3(ix,0_li,iz))= recvbuff_DOWN(4*Nz*Nx+iz*Nx+ix)
        T6(UROLL3(ix,0_li,iz))= recvbuff_DOWN(5*Nz*Nx+iz*Nx+ix)
    END DO
    END DO
    
    DO iz = 1, Nz - 2
        DO iy = 1, Ny - 2
            DO ix = 1, Nx - 2
                thiscell=UROLL3(ix, iy, iz)
                IF (w1(thiscell, x) .EQ. 1.0_dp .AND. w1(thiscell, y) .EQ. 1.0_dp .AND. w1(thiscell, z) .EQ. 1.0_dp) THEN
                
!~                 w1total=(w1(thiscell, x)+w1(thiscell, y)+w1(thiscell, z))/3
                
!~                 Vx(thiscell) = Vx(thiscell) * w1total + rho_inv *(&
                Vx(thiscell) = Vx(thiscell) + rho_inv *(&
                w2(thiscell, x) * (T1(UROLL3(ix + 1, iy, iz)) - T1(thiscell))+&
                w2(thiscell, y) * (T6(thiscell) - T6(UROLL3(ix, iy - 1, iz)))+&
                w2(thiscell, z) * (T5(thiscell) - T5(UROLL3(ix, iy, iz - 1))))

                Vy(thiscell) = Vy(thiscell) + rho_inv * (&
                w2(thiscell, x) * (T6(thiscell) - T6(UROLL3(ix - 1, iy, iz)))+&
                w2(thiscell, y) * (T2(UROLL3(ix, iy + 1, iz)) - T2(thiscell))+&
                w2(thiscell, z) * (T4(thiscell) - T4(UROLL3(ix, iy, iz - 1))))

                Vz(thiscell) = Vz(thiscell) + rho_inv * (&
                w2(thiscell, x) * (T5(thiscell) - T5(UROLL3(ix - 1, iy, iz)))+&
                w2(thiscell, y) * (T4(thiscell) - T4(UROLL3(ix, iy - 1, iz)))+&
                w2(thiscell, z) * (T3(UROLL3(ix, iy, iz + 1)) - T3(thiscell)))
                
                ELSE
                Vx_x(thiscell) = Vx_x(thiscell) * w1(thiscell, x) + &
                rho_inv * w2(thiscell, x) * (T1(UROLL3(ix + 1, iy, iz)) - T1(thiscell))

                Vx_y(thiscell) = Vx_y(thiscell) * w1(thiscell, y) + &
                rho_inv * w2(thiscell, y) * (T6(thiscell) - T6(UROLL3(ix, iy - 1, iz)))

                Vx_z(thiscell) = Vx_z(thiscell) * w1(thiscell, z) + &
                rho_inv * w2(thiscell, z) * (T5(thiscell) - T5(UROLL3(ix, iy, iz - 1)))

                Vx(thiscell) = Vx_x(thiscell) + Vx_y(thiscell) + Vx_z(thiscell)


                Vy_x(thiscell) = Vy_x(thiscell) * w1(thiscell, x) + &
                rho_inv * w2(thiscell, x) * (T6(thiscell) - T6(UROLL3(ix - 1, iy, iz)))

                Vy_y(thiscell) = Vy_y(thiscell) * w1(thiscell, y) + &
                rho_inv * w2(thiscell, y) * (T2(UROLL3(ix, iy + 1, iz)) - T2(thiscell))

                Vy_z(thiscell) = Vy_z(thiscell) * w1(thiscell, z) + &
                rho_inv * w2(thiscell, z) * (T4(thiscell) - T4(UROLL3(ix, iy, iz - 1)))

                Vy(thiscell) = Vy_x(thiscell) + Vy_y(thiscell) + Vy_z(thiscell)


                Vz_x(thiscell) = Vz_x(thiscell) * w1(thiscell, x) + &
                rho_inv * w2(thiscell, x) * (T5(thiscell) - T5(UROLL3(ix - 1, iy, iz)))

                Vz_y(thiscell) = Vz_y(thiscell) * w1(thiscell, y) + &
                rho_inv * w2(thiscell, y) * (T4(thiscell) - T4(UROLL3(ix, iy - 1, iz)))

                Vz_z(thiscell) = Vz_z(thiscell) * w1(thiscell, z) + &
                rho_inv * w2(thiscell, z) * (T3(UROLL3(ix, iy, iz + 1)) - T3(thiscell))

                Vz(thiscell) = Vz_x(thiscell) + Vz_y(thiscell) + Vz_z(thiscell)
                END IF

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

SUBROUTINE dot_source()
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE

    INTEGER(Long) :: UROLL3

    INTEGER :: ix, iy, iz

	IF (procsx.EQ.(Nprocsx-1)/2 .AND. procsy.EQ.(Nprocsy-1)/2) THEN
	    ix = (Nx/2-2)*MOD(Nprocsx-1,2)+Nx/2
	    iy = (Ny/2-2)*MOD(Nprocsy-1,2)+Ny/2
	    iz = Nz/2
!~ 	write(*,*) me, ix, iy, iz
!~ 		DO iy=1, Ny-2
		    Vz(UROLL3(ix, iy, iz)) = Vz(UROLL3(ix, iy, iz)) + &
		    (step * dt - 3 * PWIDTH)/(3 * PWIDTH) * exp(-1.0 * ((step * dt - 3.0 * PWIDTH)/(PWIDTH))**2)
!~ 	    END DO
	
	    !WRITE(*, *) (step * dt - 3 * PWIDTH)/(3 * PWIDTH) * exp(-1.0 * ((step * dt - 3.0 * PWIDTH)/(PWIDTH))**2)
	END IF
	
END SUBROUTINE dot_source

SUBROUTINE T_half_step()
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE

    INTEGER(Long) :: UROLL3
    INTEGER(Long) :: thiscell, anothercell

    INTEGER :: ix, iy, iz
    
    REAL (Double) :: w1total

    REAL(Double) :: dVxdx = 0.0_dp
    REAL(Double) :: dVydy = 0.0_dp
    REAL(Double) :: dVzdz = 0.0_dp
    
    REAL(Double) :: dVzdy = 0.0_dp
    REAL(Double) :: dVydz = 0.0_dp
    
    REAL(Double) :: dVzdx = 0.0_dp
    REAL(Double) :: dVxdz = 0.0_dp
    
    REAL(Double) :: dVydx = 0.0_dp
    REAL(Double) :: dVxdy = 0.0_dp
    
    
    DO iy=0, Ny-1
	DO ix=0, Nx-1
		iz=0
		thiscell=UROLL3(ix, iy, iz)
		iz=Nz-2
		anothercell=UROLL3(ix, iy, iz)
		Vx(thiscell)=Vx(anothercell)
		Vy(thiscell)=Vy(anothercell)
		Vz(thiscell)=Vz(anothercell)
		iz=Nz-1
		thiscell=UROLL3(ix, iy, iz)
		iz=1
		anothercell=UROLL3(ix, iy, iz)
		Vx(thiscell)=Vx(anothercell)
		Vy(thiscell)=Vy(anothercell)
		Vz(thiscell)=Vz(anothercell)
    END DO
    END DO
    
    CALL MPI_WAITALL(4, reqs, stats, ierr)
    
    DO iz=0, Nz-1
    DO iy=0, Ny-1
        Vx(UROLL3(Nx-1,iy,iz))=recvbuff_RIGHT(0*Nz*Ny+iz*Ny+iy)
        Vy(UROLL3(Nx-1,iy,iz))=recvbuff_RIGHT(1*Nz*Ny+iz*Ny+iy)
        Vz(UROLL3(Nx-1,iy,iz))=recvbuff_RIGHT(2*Nz*Ny+iz*Ny+iy)
        Vx(UROLL3(0_li,iy,iz))= recvbuff_LEFT(0*Nz*Ny+iz*Ny+iy)
        Vy(UROLL3(0_li,iy,iz))= recvbuff_LEFT(1*Nz*Ny+iz*Ny+iy)
        Vz(UROLL3(0_li,iy,iz))= recvbuff_LEFT(2*Nz*Ny+iz*Ny+iy)
    END DO
    END DO
    
    DO iz=0, Nz-1
    DO ix=0, Nx-1
        Vx(UROLL3(ix,Ny-1,iz))=recvbuff_UP(0*Nz*Nx+iz*Nx+ix)
        Vy(UROLL3(ix,Ny-1,iz))=recvbuff_UP(1*Nz*Nx+iz*Nx+ix)
        Vz(UROLL3(ix,Ny-1,iz))=recvbuff_UP(2*Nz*Nx+iz*Nx+ix)
        Vx(UROLL3(ix,0_li,iz))= recvbuff_DOWN(0*Nz*Nx+iz*Nx+ix)
        Vy(UROLL3(ix,0_li,iz))= recvbuff_DOWN(1*Nz*Nx+iz*Nx+ix)
        Vz(UROLL3(ix,0_li,iz))= recvbuff_DOWN(2*Nz*Nx+iz*Nx+ix)
    END DO
    END DO

    DO iz = 1, Nz - 2
        DO iy = 1, Ny - 2
            DO ix = 1, Nx - 2
                thiscell=UROLL3(ix, iy, iz)
                w1total=(w1(thiscell, x)+w1(thiscell, y)+w1(thiscell, z))/3
                

                dVxdx = Vx(thiscell) - Vx(UROLL3(ix - 1, iy, iz))
                dVydy = Vy(thiscell) - Vy(UROLL3(ix, iy - 1, iz))
                dVzdz = Vz(thiscell) - Vz(UROLL3(ix, iy, iz - 1))
                
                dVzdy = (Vz(UROLL3(ix, iy + 1, iz)) - Vz(thiscell))
                dVydz = (Vy(UROLL3(ix, iy, iz + 1)) - Vy(thiscell))
                
                dVzdx = (Vz(UROLL3(ix + 1, iy, iz)) - Vz(thiscell))
                dVxdz = (Vx(UROLL3(ix, iy, iz + 1)) - Vx(thiscell))
                
                dVydx = (Vy(UROLL3(ix + 1, iy, iz)) - Vy(thiscell))
                dVxdy = (Vx(UROLL3(ix, iy + 1, iz)) - Vx(thiscell))

                IF (w1(thiscell, x) .EQ. 1.0_dp .AND. w1(thiscell, y) .EQ. 1.0_dp .AND. w1(thiscell, z) .EQ. 1.0_dp) THEN
                
                T1(thiscell) = T1(thiscell) +&
                w2(thiscell, x) * c_E(1)*(dVxdx)+&
                w2(thiscell, y) * c_E(2)*(dVydy)+&
                w2(thiscell, z) * c_E(3)*(dVzdz)


                T2(thiscell) = T2(thiscell) +&
                w2(thiscell, x) * c_E(6 + 1)*(dVxdx)+&
                w2(thiscell, y) * c_E(6 + 2)*(dVydy)+&
                w2(thiscell, z) * c_E(6 + 3)*(dVzdz)

                T3(thiscell) = T3(thiscell) +&
                w2(thiscell, x) * c_E(12 + 1)*(dVxdx)+&
                w2(thiscell, y) * c_E(12 + 2)*(dVydy)+&
                w2(thiscell, z) * c_E(12 + 3)*(dVzdz)
                
                phase=(step*dt-3*PWIDTH)/(3*PWIDTH)*exp(-1.0*((step*dt-3.0*PWIDTH)/(PWIDTH))**2)
                
                dDx(thiscell)=D0x(thiscell)*phase
                dDy(thiscell)=D0y(thiscell)*phase
                dDz(thiscell)=D0z(thiscell)*phase

				dEx(thiscell)=(dDx(thiscell)-(dVzdy/deltay + dVydz/deltaz)*e_piezo(10))*beta_s(1)
				dEy(thiscell)=(dDy(thiscell)-(dVzdx/deltax + dVxdz/deltaz)*e_piezo(14))*beta_s(5)
				dEz(thiscell)=(dDz(thiscell)-(dVydx/deltax + dVxdy/deltay)*e_piezo(18))*beta_s(9)
				
                T4(thiscell) = T4(thiscell) +&
                w2(thiscell, y) * c_E(18 + 4)*dVzdy+&
                w2(thiscell, z) * c_E(18 + 4)*dVydz+ &
                dEx(thiscell)*e_piezo(10)*dt

                T5(thiscell) = T5(thiscell) +&
                w2(thiscell, x) * c_E(24 + 5)*dVzdx+&
                w2(thiscell, z) * c_E(24 + 5)*dVxdz+ &
                dEy(thiscell)*e_piezo(14)*dt

                T6(thiscell) = T6(thiscell) +&
                w2(thiscell, x) * c_E(30 + 6)*dVydx+&
                w2(thiscell, y) * c_E(30 + 6)*dVxdy+ &
                dEz(thiscell)*e_piezo(18)*dt
                
                ELSE
                
                T1_x(thiscell) = T1_x(thiscell) * w1(thiscell, x) +&
                w2(thiscell, x) * c_E(1)*(dVxdx)
                T1_y(thiscell) = T1_y(thiscell) * w1(thiscell, y) +&
                w2(thiscell, y) * c_E(2)*(dVydy)
                T1_z(thiscell) = T1_z(thiscell) * w1(thiscell, z) +&
                w2(thiscell, z) * c_E(3)*(dVzdz)
                T1(thiscell) = T1_x(thiscell) + T1_y(thiscell) + T1_z(thiscell)


                T2_x(thiscell) = T2_x(thiscell) * w1(thiscell, x) +&
                w2(thiscell, x) * c_E(6 + 1)*(dVxdx)
                T2_y(thiscell) = T2_y(thiscell) * w1(thiscell, y) +&
                w2(thiscell, y) * c_E(6 + 2)*(dVydy)
                T2_z(thiscell) = T2_z(thiscell) * w1(thiscell, z) +&
                w2(thiscell, z) * c_E(6 + 3)*(dVzdz)
                T2(thiscell) = T2_x(thiscell) + T2_y(thiscell) + T2_z(thiscell)

                T3_x(thiscell) = T3_x(thiscell) * w1(thiscell, x) +&
                w2(thiscell, x) * c_E(12 + 1)*(dVxdx)
                T3_y(thiscell) = T3_y(thiscell) * w1(thiscell, y) +&
                w2(thiscell, y) * c_E(12 + 2)*(dVydy)
                T3_z(thiscell) = T3_z(thiscell) * w1(thiscell, z) +&
                w2(thiscell, z) * c_E(12 + 3)*(dVzdz)
                T3(thiscell) = T3_x(thiscell) + T3_y(thiscell) + T3_z(thiscell)
                
                phase=(step*dt-3*PWIDTH)/(3*PWIDTH)*exp(-1.0*((step*dt-3.0*PWIDTH)/(PWIDTH))**2)
                
                dDx(thiscell)=D0x(thiscell)*phase
                dDy(thiscell)=D0y(thiscell)*phase
                dDz(thiscell)=D0z(thiscell)*phase

				dEx(thiscell)=(dDx(thiscell)-(dVzdy/deltay + dVydz/deltaz)*e_piezo(10))*beta_s(1)
				dEy(thiscell)=(dDy(thiscell)-(dVzdx/deltax + dVxdz/deltaz)*e_piezo(14))*beta_s(5)
				dEz(thiscell)=(dDz(thiscell)-(dVydx/deltax + dVxdy/deltay)*e_piezo(18))*beta_s(9)
				
                T4_y(thiscell) = T4_y(thiscell) * w1(thiscell, y) +&
                w2(thiscell, y) * c_E(18 + 4)*dVzdy
                T4_z(thiscell) = T4_z(thiscell) * w1(thiscell, z) +&
                w2(thiscell, z) * c_E(18 + 4)*dVydz
                T4(thiscell) = T4_y(thiscell) + T4_z(thiscell) + &
                dEx(thiscell)*e_piezo(10)*dt

                T5_x(thiscell) = T5_x(thiscell) * w1(thiscell, x) +&
                w2(thiscell, x) * c_E(24 + 5)*dVzdx
                T5_z(thiscell) = T5_z(thiscell) * w1(thiscell, z) +&
                w2(thiscell, z) * c_E(24 + 5)*dVxdz
                T5(thiscell) = T5_x(thiscell) + T5_z(thiscell) + &
                dEy(thiscell)*e_piezo(14)*dt

                T6_x(thiscell) = T6_x(thiscell) * w1(thiscell, x) +&
                w2(thiscell, x) * c_E(30 + 6)*dVydx
                T6_y(thiscell) = T6_y(thiscell) * w1(thiscell, y) +&
                w2(thiscell, y) * c_E(30 + 6)*dVxdy
                T6(thiscell) = T6_x(thiscell) + T6_y(thiscell) + &
                dEz(thiscell)*e_piezo(18)*dt
                END IF

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


SUBROUTINE share_v()
	USE mpi
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
	IMPLICIT NONE
	INTEGER(Long) :: UROLL3
	INTEGER(Short) UROLLPROC
	INTEGER :: ix, iy, iz, nextprocid
    
    CALL MPI_IRECV(recvbuff_RIGHT, vbuffsizex, MPI_DOUBLE_PRECISION, &
                           nRIGHT, 1, MPI_COMM_WORLD, reqs(1), ierr)
    CALL MPI_IRECV(recvbuff_LEFT , vbuffsizex, MPI_DOUBLE_PRECISION, &
                           nLEFT , 3, MPI_COMM_WORLD, reqs(2), ierr)
                           
    CALL MPI_IRECV(recvbuff_DOWN , vbuffsizey, MPI_DOUBLE_PRECISION, &
                           nDOWN , 2, MPI_COMM_WORLD, reqs(3), ierr)
    CALL MPI_IRECV(recvbuff_UP   , vbuffsizey, MPI_DOUBLE_PRECISION, &
                           nUP   , 0, MPI_COMM_WORLD, reqs(4), ierr)
                   
    DO iz=0, Nz-1
    DO iy=0, Ny-1
        sendbuff_RIGHT(0*Nz*Ny+iz*Ny+iy)= Vx(UROLL3(Nx-2,iy,iz))
        sendbuff_RIGHT(1*Nz*Ny+iz*Ny+iy)= Vy(UROLL3(Nx-2,iy,iz))
        sendbuff_RIGHT(2*Nz*Ny+iz*Ny+iy)= Vz(UROLL3(Nx-2,iy,iz))
        sendbuff_LEFT(0*Nz*Ny+iz*Ny+iy) = Vx(UROLL3(1_li,iy,iz))
        sendbuff_LEFT(1*Nz*Ny+iz*Ny+iy) = Vy(UROLL3(1_li,iy,iz))
        sendbuff_LEFT(2*Nz*Ny+iz*Ny+iy) = Vz(UROLL3(1_li,iy,iz))
    END DO
    END DO
        
    DO iz=0, Nz-1
    DO ix=0, Nx-1
        sendbuff_UP(0*Nz*Nx+iz*Nx+ix)= Vx(UROLL3(ix,Ny-2,iz))
        sendbuff_UP(1*Nz*Nx+iz*Nx+ix)= Vy(UROLL3(ix,Ny-2,iz))
        sendbuff_UP(2*Nz*Nx+iz*Nx+ix)= Vz(UROLL3(ix,Ny-2,iz))
        sendbuff_DOWN(0*Nz*Nx+iz*Nx+ix) = Vx(UROLL3(ix,1_li,iz))
        sendbuff_DOWN(1*Nz*Nx+iz*Nx+ix) = Vy(UROLL3(ix,1_li,iz))
        sendbuff_DOWN(2*Nz*Nx+iz*Nx+ix) = Vz(UROLL3(ix,1_li,iz))
    END DO
    END DO
    
    CALL MPI_ISEND(sendbuff_LEFT , vbuffsizex, MPI_DOUBLE_PRECISION, &
                           nLEFT , 1, MPI_COMM_WORLD, reqs(5), ierr)
    CALL MPI_ISEND(sendbuff_RIGHT, vbuffsizex, MPI_DOUBLE_PRECISION, &
                           nRIGHT, 3, MPI_COMM_WORLD, reqs(6), ierr)
                           
    CALL MPI_ISEND(sendbuff_DOWN , vbuffsizey, MPI_DOUBLE_PRECISION, &
                           nDOWN , 0, MPI_COMM_WORLD, reqs(7), ierr)
    CALL MPI_ISEND(sendbuff_UP   , vbuffsizey, MPI_DOUBLE_PRECISION, &
                           nUP   , 2, MPI_COMM_WORLD, reqs(8), ierr)

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
	
	
    CALL MPI_IRECV(recvbuff_RIGHT, 2*vbuffsizex, MPI_DOUBLE_PRECISION, &
                           nRIGHT, 1, MPI_COMM_WORLD, reqs(1), ierr)
    CALL MPI_IRECV(recvbuff_LEFT , 2*vbuffsizex, MPI_DOUBLE_PRECISION, &
                           nLEFT , 3, MPI_COMM_WORLD, reqs(2), ierr)
                           
    CALL MPI_IRECV(recvbuff_DOWN , 2*vbuffsizey, MPI_DOUBLE_PRECISION, &
                           nDOWN , 2, MPI_COMM_WORLD, reqs(3), ierr)
    CALL MPI_IRECV(recvbuff_UP   , 2*vbuffsizey, MPI_DOUBLE_PRECISION, &
                           nUP   , 0, MPI_COMM_WORLD, reqs(4), ierr)
                   
    DO iz=0, Nz-1
    DO iy=0, Ny-1
        sendbuff_RIGHT(0*Nz*Ny+iz*Ny+iy)= T1(UROLL3(Nx-2,iy,iz))
        sendbuff_RIGHT(1*Nz*Ny+iz*Ny+iy)= T2(UROLL3(Nx-2,iy,iz))
        sendbuff_RIGHT(2*Nz*Ny+iz*Ny+iy)= T3(UROLL3(Nx-2,iy,iz))
        sendbuff_RIGHT(3*Nz*Ny+iz*Ny+iy)= T4(UROLL3(Nx-2,iy,iz))
        sendbuff_RIGHT(4*Nz*Ny+iz*Ny+iy)= T5(UROLL3(Nx-2,iy,iz))
        sendbuff_RIGHT(5*Nz*Ny+iz*Ny+iy)= T6(UROLL3(Nx-2,iy,iz))
        sendbuff_LEFT(0*Nz*Ny+iz*Ny+iy) = T1(UROLL3(1_li,iy,iz))
        sendbuff_LEFT(1*Nz*Ny+iz*Ny+iy) = T2(UROLL3(1_li,iy,iz))
        sendbuff_LEFT(2*Nz*Ny+iz*Ny+iy) = T3(UROLL3(1_li,iy,iz))
        sendbuff_LEFT(3*Nz*Ny+iz*Ny+iy) = T4(UROLL3(1_li,iy,iz))
        sendbuff_LEFT(4*Nz*Ny+iz*Ny+iy) = T5(UROLL3(1_li,iy,iz))
        sendbuff_LEFT(5*Nz*Ny+iz*Ny+iy) = T6(UROLL3(1_li,iy,iz))
    END DO
    END DO
        
    DO iz=0, Nz-1
    DO ix=0, Nx-1
        sendbuff_UP(0*Nz*Nx+iz*Nx+ix)= T1(UROLL3(ix,Ny-2,iz))
        sendbuff_UP(1*Nz*Nx+iz*Nx+ix)= T2(UROLL3(ix,Ny-2,iz))
        sendbuff_UP(2*Nz*Nx+iz*Nx+ix)= T3(UROLL3(ix,Ny-2,iz))
        sendbuff_UP(3*Nz*Nx+iz*Nx+ix)= T4(UROLL3(ix,Ny-2,iz))
        sendbuff_UP(4*Nz*Nx+iz*Nx+ix)= T5(UROLL3(ix,Ny-2,iz))
        sendbuff_UP(5*Nz*Nx+iz*Nx+ix)= T6(UROLL3(ix,Ny-2,iz))
        sendbuff_DOWN(0*Nz*Nx+iz*Nx+ix) = T1(UROLL3(ix,1_li,iz))
        sendbuff_DOWN(1*Nz*Nx+iz*Nx+ix) = T2(UROLL3(ix,1_li,iz))
        sendbuff_DOWN(2*Nz*Nx+iz*Nx+ix) = T3(UROLL3(ix,1_li,iz))
        sendbuff_DOWN(3*Nz*Nx+iz*Nx+ix) = T4(UROLL3(ix,1_li,iz))
        sendbuff_DOWN(4*Nz*Nx+iz*Nx+ix) = T5(UROLL3(ix,1_li,iz))
        sendbuff_DOWN(5*Nz*Nx+iz*Nx+ix) = T6(UROLL3(ix,1_li,iz))
    END DO
    END DO
    
    CALL MPI_ISEND(sendbuff_LEFT , 2*vbuffsizex, MPI_DOUBLE_PRECISION, &
                           nLEFT , 1, MPI_COMM_WORLD, reqs(5), ierr)
    CALL MPI_ISEND(sendbuff_RIGHT, 2*vbuffsizex, MPI_DOUBLE_PRECISION, &
                           nRIGHT, 3, MPI_COMM_WORLD, reqs(6), ierr)
                           
    CALL MPI_ISEND(sendbuff_DOWN , 2*vbuffsizey, MPI_DOUBLE_PRECISION, &
                           nDOWN , 0, MPI_COMM_WORLD, reqs(7), ierr)
    CALL MPI_ISEND(sendbuff_UP   , 2*vbuffsizey, MPI_DOUBLE_PRECISION, &
                           nUP   , 2, MPI_COMM_WORLD, reqs(8), ierr)
	    
END SUBROUTINE share_T





SUBROUTINE Kinetic_Energy()
	USE mpi
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
	IMPLICIT NONE
	INTEGER(Long) :: UROLL3
	INTEGER(Short) UROLLPROC
	INTEGER :: ix, iy, iz, nextprocid
    INTEGER(Long) :: thiscell
    
    U_k=0;
    DO iz = 1, Nz - 2
        DO iy = 1, Ny - 2
            DO ix = 1, Nx - 2
                thiscell=UROLL3(ix, iy, iz)
                U_k=U_k+Vx(thiscell)**2+Vy(thiscell)**2+Vz(thiscell)**2
            END DO
        END DO
    END DO
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(U_k, U_k_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    
END SUBROUTINE Kinetic_Energy





