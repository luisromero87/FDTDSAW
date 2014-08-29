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
    INTEGER(Long) :: thiscell


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
    
    REAL(Double) :: dVzdy = 0.0_dp
    REAL(Double) :: dVydz = 0.0_dp
    
    REAL(Double) :: dVzdx = 0.0_dp
    REAL(Double) :: dVxdz = 0.0_dp
    
    REAL(Double) :: dVydx = 0.0_dp
    REAL(Double) :: dVxdy = 0.0_dp
    

    DO iz = 1, Nz - 2
        DO iy = 1, Ny - 2
            DO ix = 1, Nx - 2

                dVxdx = Vx(UROLL3(ix, iy, iz)) - Vx(UROLL3(ix - 1, iy, iz))
                dVydy = Vy(UROLL3(ix, iy, iz)) - Vy(UROLL3(ix, iy - 1, iz))
                dVzdz = Vz(UROLL3(ix, iy, iz)) - Vz(UROLL3(ix, iy, iz - 1))
                
                dVzdy = (Vz(UROLL3(ix, iy + 1, iz)) - Vz(UROLL3(ix, iy, iz)))
                dVydz = (Vy(UROLL3(ix, iy, iz + 1)) - Vy(UROLL3(ix, iy, iz)))
                
                dVzdx = (Vz(UROLL3(ix + 1, iy, iz)) - Vz(UROLL3(ix, iy, iz)))
                dVxdz = (Vx(UROLL3(ix, iy, iz + 1)) - Vx(UROLL3(ix, iy, iz)))
                
                dVydx = (Vy(UROLL3(ix + 1, iy, iz)) - Vy(UROLL3(ix, iy, iz)))
                dVxdy = (Vx(UROLL3(ix, iy + 1, iz)) - Vx(UROLL3(ix, iy, iz)))

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
                
                dDx(UROLL3(ix, iy, iz))=D0x(UROLL3(ix, iy, iz))*phase
                dDy(UROLL3(ix, iy, iz))=D0y(UROLL3(ix, iy, iz))*phase
                dDz(UROLL3(ix, iy, iz))=D0z(UROLL3(ix, iy, iz))*phase

				dEx(UROLL3(ix, iy, iz))=(dDx(UROLL3(ix, iy, iz))-(dVzdy/deltay + dVydz/deltaz)*e_piezo(10))*beta_s(1)
				dEy(UROLL3(ix, iy, iz))=(dDy(UROLL3(ix, iy, iz))-(dVzdx/deltax + dVxdz/deltaz)*e_piezo(14))*beta_s(5)
				dEz(UROLL3(ix, iy, iz))=(dDz(UROLL3(ix, iy, iz))-(dVydx/deltax + dVxdy/deltay)*e_piezo(18))*beta_s(9)
				
                T4_y(UROLL3(ix, iy, iz)) = T4_y(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), y) +&
                w2(UROLL3(ix, iy, iz), y) * c_E(18 + 4)*dVzdy
                T4_z(UROLL3(ix, iy, iz)) = T4_z(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), z) +&
                w2(UROLL3(ix, iy, iz), z) * c_E(18 + 4)*dVydz
                T4(UROLL3(ix, iy, iz)) = T4_y(UROLL3(ix, iy, iz)) + T4_z(UROLL3(ix, iy, iz)) + &
                dEx(UROLL3(ix, iy, iz))*e_piezo(10)*dt

                T5_x(UROLL3(ix, iy, iz)) = T5_x(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), x) +&
                w2(UROLL3(ix, iy, iz), x) * c_E(24 + 5)*dVzdx
                T5_z(UROLL3(ix, iy, iz)) = T5_z(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), z) +&
                w2(UROLL3(ix, iy, iz), z) * c_E(24 + 5)*dVxdz
                T5(UROLL3(ix, iy, iz)) = T5_x(UROLL3(ix, iy, iz)) + T5_z(UROLL3(ix, iy, iz)) + &
                dEy(UROLL3(ix, iy, iz))*e_piezo(14)*dt

                T6_x(UROLL3(ix, iy, iz)) = T6_x(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), x) +&
                w2(UROLL3(ix, iy, iz), x) * c_E(30 + 6)*dVydx
                T6_y(UROLL3(ix, iy, iz)) = T6_y(UROLL3(ix, iy, iz)) * w1(UROLL3(ix, iy, iz), y) +&
                w2(UROLL3(ix, iy, iz), y) * c_E(30 + 6)*dVxdy
                T6(UROLL3(ix, iy, iz)) = T6_x(UROLL3(ix, iy, iz)) + T6_y(UROLL3(ix, iy, iz)) + &
                dEz(UROLL3(ix, iy, iz))*e_piezo(18)*dt

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

	IF (procsx.EQ.(Nprocsx-1)/2 .AND. procsy.EQ.(Nprocsy-1)/2) THEN
	    ix = (Nx/2-2)*MOD(Nprocsx-1,2)+Nx/2
	    iy = (Ny/2-2)*MOD(Nprocsy-1,2)+Ny/2
	    iz = 1
!~ 	write(*,*) me, ix, iy, iz
!~ 		DO iy=1, Ny-2
		    Vz(UROLL3(ix, iy, iz)) = Vz(UROLL3(ix, iy, iz)) + &
		    (step * dt - 3 * PWIDTH)/(3 * PWIDTH) * exp(-1.0 * ((step * dt - 3.0 * PWIDTH)/(PWIDTH))**2)
!~ 	    END DO
	
	    !WRITE(*, *) (step * dt - 3 * PWIDTH)/(3 * PWIDTH) * exp(-1.0 * ((step * dt - 3.0 * PWIDTH)/(PWIDTH))**2)
	END IF
	
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
	    
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
		

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
	    
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

END SUBROUTINE share_T









