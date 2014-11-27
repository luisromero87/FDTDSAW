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
		T1(ix, iy, 0)=T1(ix, iy, Nz-2)
		T2(ix, iy, 0)=T2(ix, iy, Nz-2)
		T3(ix, iy, 0)=T3(ix, iy, Nz-2)
		T4(ix, iy, 0)=T4(ix, iy, Nz-2)
		T5(ix, iy, 0)=T5(ix, iy, Nz-2)
		T6(ix, iy, 0)=T6(ix, iy, Nz-2)
		iz=Nz-1
		thiscell=UROLL3(ix, iy, iz)
		iz=1
		anothercell=UROLL3(ix, iy, iz)
		T1(ix, iy, Nz-1)=T1(ix, iy, 1)
		T2(ix, iy, Nz-1)=T2(ix, iy, 1)
		T3(ix, iy, Nz-1)=T3(ix, iy, 1)
		T4(ix, iy, Nz-1)=T4(ix, iy, 1)
		T5(ix, iy, Nz-1)=T5(ix, iy, 1)
		T6(ix, iy, Nz-1)=T6(ix, iy, 1)
    END DO
    END DO
        
    CALL MPI_WAITALL(4, reqs, stats, ierr)
    
    DO iz=0, Nz-1
    DO iy=0, Ny-1
        T1(Nx-1,iy,iz)=recvbuff_RIGHT(0*Nz*Ny+iz*Ny+iy)
        T2(Nx-1,iy,iz)=recvbuff_RIGHT(1*Nz*Ny+iz*Ny+iy)
        T3(Nx-1,iy,iz)=recvbuff_RIGHT(2*Nz*Ny+iz*Ny+iy)
        T4(Nx-1,iy,iz)=recvbuff_RIGHT(3*Nz*Ny+iz*Ny+iy)
        T5(Nx-1,iy,iz)=recvbuff_RIGHT(4*Nz*Ny+iz*Ny+iy)
        T6(Nx-1,iy,iz)=recvbuff_RIGHT(5*Nz*Ny+iz*Ny+iy)
        T1(0_li,iy,iz)= recvbuff_LEFT(0*Nz*Ny+iz*Ny+iy)
        T2(0_li,iy,iz)= recvbuff_LEFT(1*Nz*Ny+iz*Ny+iy)
        T3(0_li,iy,iz)= recvbuff_LEFT(2*Nz*Ny+iz*Ny+iy)
        T4(0_li,iy,iz)= recvbuff_LEFT(3*Nz*Ny+iz*Ny+iy)
        T5(0_li,iy,iz)= recvbuff_LEFT(4*Nz*Ny+iz*Ny+iy)
        T6(0_li,iy,iz)= recvbuff_LEFT(5*Nz*Ny+iz*Ny+iy)
    END DO
    END DO
    
    DO iz=0, Nz-1
    DO ix=0, Nx-1
        T1(ix,Ny-1,iz)=recvbuff_UP(0*Nz*Nx+iz*Nx+ix)
        T2(ix,Ny-1,iz)=recvbuff_UP(1*Nz*Nx+iz*Nx+ix)
        T3(ix,Ny-1,iz)=recvbuff_UP(2*Nz*Nx+iz*Nx+ix)
        T4(ix,Ny-1,iz)=recvbuff_UP(3*Nz*Nx+iz*Nx+ix)
        T5(ix,Ny-1,iz)=recvbuff_UP(4*Nz*Nx+iz*Nx+ix)
        T6(ix,Ny-1,iz)=recvbuff_UP(5*Nz*Nx+iz*Nx+ix)
        T1(ix,0_li,iz)= recvbuff_DOWN(0*Nz*Nx+iz*Nx+ix)
        T2(ix,0_li,iz)= recvbuff_DOWN(1*Nz*Nx+iz*Nx+ix)
        T3(ix,0_li,iz)= recvbuff_DOWN(2*Nz*Nx+iz*Nx+ix)
        T4(ix,0_li,iz)= recvbuff_DOWN(3*Nz*Nx+iz*Nx+ix)
        T5(ix,0_li,iz)= recvbuff_DOWN(4*Nz*Nx+iz*Nx+ix)
        T6(ix,0_li,iz)= recvbuff_DOWN(5*Nz*Nx+iz*Nx+ix)
    END DO
    END DO
    
    DO iz = 1, Nz - 2
        DO iy = 1, Ny - 2
            DO ix = 1, Nx - 2
                thiscell=UROLL3(ix, iy, iz)
                IF (w1(ix, iy, iz, x) .EQ. 1.0_dp .AND. w1(ix, iy, iz, y) .EQ. 1.0_dp .AND. w1(ix, iy, iz, z) .EQ. 1.0_dp) THEN
                
!~                 w1total=(w1(ix, iy, iz, x)+w1(ix, iy, iz, y)+w1(ix, iy, iz, z))/3
                
!~                 Vx(ix, iy, iz) = Vx(ix, iy, iz) * w1total + rho_inv *(&
                Vx(ix,iy,iz) = Vx(ix,iy,iz) + rho_inv *(&
                w2(ix,iy,iz,x) * (T1(ix + 1,iy,iz) - T1(ix,iy,iz))+&
                w2(ix,iy,iz,y) * (T6(ix,iy,iz) - T6(ix,iy - 1,iz))+&
                w2(ix,iy,iz,z) * (T5(ix,iy,iz) - T5(ix,iy,iz - 1)))

                Vy(ix,iy,iz) = Vy(ix,iy,iz) + rho_inv * (&
                w2(ix,iy,iz,x) * (T6(ix,iy,iz) - T6(ix - 1,iy,iz))+&
                w2(ix,iy,iz,y) * (T2(ix,iy + 1,iz) - T2(ix,iy,iz))+&
                w2(ix,iy,iz,z) * (T4(ix,iy,iz) - T4(ix,iy,iz - 1)))

                Vz(ix,iy,iz) = Vz(ix,iy,iz) + rho_inv * (&
                w2(ix,iy,iz,x) * (T5(ix,iy,iz) - T5(ix - 1,iy,iz))+&
                w2(ix,iy,iz,y) * (T4(ix,iy,iz) - T4(ix,iy - 1,iz))+&
                w2(ix,iy,iz,z) * (T3(ix,iy,iz + 1) - T3(ix,iy,iz)))
                
                ELSE
                Vx_x(ix,iy,iz) = Vx_x(ix,iy,iz) * w1(ix,iy,iz,x) + rho_inv * &
                w2(ix,iy,iz,x) * (T1(ix + 1,iy,iz) - T1(ix,iy,iz))

                Vx_y(ix,iy,iz) = Vx_y(ix,iy,iz) * w1(ix,iy,iz,y) + rho_inv * &
                w2(ix,iy,iz,y) * (T6(ix,iy,iz) - T6(ix,iy - 1,iz))

                Vx_z(ix,iy,iz) = Vx_z(ix,iy,iz) * w1(ix,iy,iz,z) + rho_inv * &
                w2(ix,iy,iz,z) * (T5(ix,iy,iz) - T5(ix,iy,iz - 1))

                Vx(ix,iy,iz) = Vx_x(ix,iy,iz) + Vx_y(ix,iy,iz) + Vx_z(ix,iy,iz)


                Vy_x(ix,iy,iz) = Vy_x(ix,iy,iz) * w1(ix,iy,iz,x) + rho_inv * &
                w2(ix,iy,iz,x) * (T6(ix,iy,iz) - T6(ix - 1,iy,iz))

                Vy_y(ix,iy,iz) = Vy_y(ix,iy,iz) * w1(ix,iy,iz,y) + rho_inv * &
                w2(ix,iy,iz,y) * (T2(ix,iy + 1,iz) - T2(ix,iy,iz))

                Vy_z(ix,iy,iz) = Vy_z(ix,iy,iz) * w1(ix,iy,iz,z) + rho_inv * &
                w2(ix,iy,iz,z) * (T4(ix,iy,iz) - T4(ix,iy,iz - 1))

                Vy(ix,iy,iz) = Vy_x(ix,iy,iz) + Vy_y(ix,iy,iz) + Vy_z(ix,iy,iz)


                Vz_x(ix,iy,iz) = Vz_x(ix,iy,iz) * w1(ix,iy,iz,x) + rho_inv * &
                w2(ix,iy,iz,x) * (T5(ix,iy,iz) - T5(ix - 1,iy,iz))

                Vz_y(ix,iy,iz) = Vz_y(ix,iy,iz) * w1(ix,iy,iz,y) + rho_inv * &
                w2(ix,iy,iz,y) * (T4(ix,iy,iz) - T4(ix,iy - 1,iz))

                Vz_z(ix,iy,iz) = Vz_z(ix,iy,iz) * w1(ix,iy,iz,z) + rho_inv * &
                w2(ix,iy,iz,z) * (T3(ix,iy,iz + 1) - T3(ix,iy,iz))

                Vz(ix,iy,iz) = Vz_x(ix,iy,iz) + Vz_y(ix,iy,iz) + Vz_z(ix,iy,iz)
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
            Vz(ix, iy, iz - 1) = Vz(ix, iy, iz)+(dz(iz)/c_E(3,3))*(&
            c_E(3,1)/dx(ix) * (Vx(ix, iy, iz) - Vx(ix - 1, iy, iz))+&
            c_E(3,2)/dy(iy) * (Vy(ix, iy, iz) - Vy(ix, iy - 1, iz)))
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
! 		DO iz=1, Nz-2
		    Vz(ix, iy, iz) = Vz(ix, iy, iz) + &
		    (step * dt - 3 * PWIDTH)/(3 * PWIDTH) * exp(-1.0 * ((step * dt - 3.0 * PWIDTH)/(PWIDTH))**2)
		    Vy(ix, iy, iz) = Vy(ix, iy, iz) + &
		    (step * dt - 3 * PWIDTH)/(3 * PWIDTH) * exp(-1.0 * ((step * dt - 3.0 * PWIDTH)/(PWIDTH))**2)
		    Vx(ix, iy, iz) = Vx(ix, iy, iz) + &
		    (step * dt - 3 * PWIDTH)/(3 * PWIDTH) * exp(-1.0 * ((step * dt - 3.0 * PWIDTH)/(PWIDTH))**2)
! 	    END DO
	
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
		Vx(ix, iy, 0)=Vx(ix, iy, Nz-2)
		Vy(ix, iy, 0)=Vy(ix, iy, Nz-2)
		Vz(ix, iy, 0)=Vz(ix, iy, Nz-2)
		iz=Nz-1
		thiscell=UROLL3(ix, iy, iz)
		iz=1
		anothercell=UROLL3(ix, iy, iz)
		Vx(ix, iy, Nz-1)=Vx(ix, iy, 1)
		Vy(ix, iy, Nz-1)=Vy(ix, iy, 1)
		Vz(ix, iy, Nz-1)=Vz(ix, iy, 1)
    END DO
    END DO
    
    CALL MPI_WAITALL(4, reqs, stats, ierr)
    
    DO iz=0, Nz-1
    DO iy=0, Ny-1
        Vx(Nx-1,iy,iz)=recvbuff_RIGHT(0*Nz*Ny+iz*Ny+iy)
        Vy(Nx-1,iy,iz)=recvbuff_RIGHT(1*Nz*Ny+iz*Ny+iy)
        Vz(Nx-1,iy,iz)=recvbuff_RIGHT(2*Nz*Ny+iz*Ny+iy)
        Vx(0_li,iy,iz)= recvbuff_LEFT(0*Nz*Ny+iz*Ny+iy)
        Vy(0_li,iy,iz)= recvbuff_LEFT(1*Nz*Ny+iz*Ny+iy)
        Vz(0_li,iy,iz)= recvbuff_LEFT(2*Nz*Ny+iz*Ny+iy)
    END DO
    END DO
    
    DO iz=0, Nz-1
    DO ix=0, Nx-1
        Vx(ix,Ny-1,iz)=recvbuff_UP(0*Nz*Nx+iz*Nx+ix)
        Vy(ix,Ny-1,iz)=recvbuff_UP(1*Nz*Nx+iz*Nx+ix)
        Vz(ix,Ny-1,iz)=recvbuff_UP(2*Nz*Nx+iz*Nx+ix)
        Vx(ix,0_li,iz)= recvbuff_DOWN(0*Nz*Nx+iz*Nx+ix)
        Vy(ix,0_li,iz)= recvbuff_DOWN(1*Nz*Nx+iz*Nx+ix)
        Vz(ix,0_li,iz)= recvbuff_DOWN(2*Nz*Nx+iz*Nx+ix)
    END DO
    END DO

    DO iz = 1, Nz - 2
        DO iy = 1, Ny - 2
            DO ix = 1, Nx - 2
                thiscell=UROLL3(ix, iy, iz)                

                dVxdx = Vx(ix, iy, iz) - Vx(ix - 1, iy, iz)
                dVydy = Vy(ix, iy, iz) - Vy(ix, iy - 1, iz)
                dVzdz = Vz(ix, iy, iz) - Vz(ix, iy, iz - 1)
                
                dVzdy = (Vz(ix, iy + 1, iz) - Vz(ix, iy, iz))
                dVydz = (Vy(ix, iy, iz + 1) - Vy(ix, iy, iz))
                
                dVzdx = (Vz(ix + 1, iy, iz) - Vz(ix, iy, iz))
                dVxdz = (Vx(ix, iy, iz + 1) - Vx(ix, iy, iz))
                
                dVydx = (Vy(ix + 1, iy, iz) - Vy(ix, iy, iz))
                dVxdy = (Vx(ix, iy + 1, iz) - Vx(ix, iy, iz))
                
                phase=(step*dt-3*PWIDTH)/(3*PWIDTH)*exp(-1.0*((step*dt-3.0*PWIDTH)/(PWIDTH))**2)
                
                dDx(ix, iy, iz)=D0x(ix, iy, iz)*phase
                dDy(ix, iy, iz)=D0y(ix, iy, iz)*phase
                dDz(ix, iy, iz)=D0z(ix, iy, iz)*phase
                
!                IF (procsx.EQ.(Nprocsx-1)/2 .AND. procsy.EQ.(Nprocsy-1)/2) THEN
!                IF (ix .EQ. (Nx/2-2)*MOD(Nprocsx-1,2)+Nx/2 .AND. iy .EQ. (Ny/2-2)*MOD(Nprocsy-1,2)+Ny/2 .AND. iz .EQ. Nz/2) THEN
!                  dDx(thiscell)=phase
!                ENDIF
!                ENDIF
                
                !Cubic                 
				dEx(ix, iy, iz)=(dDx(ix, iy, iz)-(dVzdy/deltay + dVydz/deltaz)*e_piezo(1,3))*beta_s(1,1)
				dEy(ix, iy, iz)=(dDy(ix, iy, iz)-(dVzdx/deltax + dVxdz/deltaz)*e_piezo(2,5))*beta_s(2,2)
				dEz(ix, iy, iz)=(dDz(ix, iy, iz)-(dVydx/deltax + dVxdy/deltay)*e_piezo(3,6))*beta_s(3,3)

                !Orthorhombic 2mm
!				dEx(thiscell)=(dDx(thiscell)-(dVxdz/deltaz + dVzdx/deltax)*e_piezo(1,5))*beta_s(1,1)
!				dEy(thiscell)=(dDy(thiscell)-(dVydz/deltaz + dVzdy/deltay)*e_piezo(2,4))*beta_s(2,2)
!				dEz(thiscell)=( dDz(thiscell) &
!                -(e_piezo(3,1)*dVxdx/deltax + e_piezo(3,2)*dVydy/deltay + e_piezo(3,3)*dVzdz/deltaz) )*beta_s(3,3)

                IF (w1(ix,iy,iz,x) .EQ. 1.0_dp .AND. w1(ix,iy,iz,y) .EQ. 1.0_dp .AND. w1(ix,iy,iz,z) .EQ. 1.0_dp) THEN
!                cubic 
                T1(ix, iy, iz) = T1(ix, iy, iz) +&
                w2(ix, iy, iz, x) * c_E(1,1)*(dVxdx)+&
                w2(ix, iy, iz, y) * c_E(1,2)*(dVydy)+&
                w2(ix, iy, iz, z) * c_E(1,3)*(dVzdz)

                T2(ix, iy, iz) = T2(ix, iy, iz) +&
                w2(ix, iy, iz, x) * c_E(2,1)*(dVxdx)+&
                w2(ix, iy, iz, y) * c_E(2,2)*(dVydy)+&
                w2(ix, iy, iz, z) * c_E(2,3)*(dVzdz)

                T3(ix, iy, iz) = T3(ix, iy, iz) +&
                w2(ix, iy, iz, x) * c_E(3,1)*(dVxdx)+&
                w2(ix, iy, iz, y) * c_E(3,2)*(dVydy)+&
                w2(ix, iy, iz, z) * c_E(3,3)*(dVzdz)
				
                T4(ix, iy, iz) = T4(ix, iy, iz) +&
                w2(ix, iy, iz, y) * c_E(4,4)*dVzdy+&
                w2(ix, iy, iz, z) * c_E(4,4)*dVydz - &
                dEx(ix, iy, iz)*e_piezo(1,4)*dt

                T5(ix, iy, iz) = T5(ix, iy, iz) +&
                w2(ix, iy, iz, x) * c_E(5,5)*dVzdx+&
                w2(ix, iy, iz, z) * c_E(5,5)*dVxdz - &
                dEy(ix, iy, iz)*e_piezo(2,5)*dt

                T6(ix, iy, iz) = T6(ix, iy, iz) +&
                w2(ix, iy, iz, x) * c_E(6,6)*dVydx+&
                w2(ix, iy, iz, y) * c_E(6,6)*dVxdy - &
                dEz(ix, iy, iz)*e_piezo(3,6)*dt
!                
                !Orthorhombic 2mm
!                T1(thiscell) = T1(thiscell) +&
!                w2(thiscell, x) * c_E(1,1)*(dVxdx)+&
!                w2(thiscell, y) * c_E(1,2)*(dVydy)+&
!                w2(thiscell, z) * c_E(1,3)*(dVzdz) - e_piezo(3,1)*dEz(thiscell)*dt
!
!                T2(thiscell) = T2(thiscell) +&
!                w2(thiscell, x) * c_E(2,1)*(dVxdx)+&
!                w2(thiscell, y) * c_E(2,2)*(dVydy)+&
!                w2(thiscell, z) * c_E(2,3)*(dVzdz) - e_piezo(3,2)*dEz(thiscell)*dt
!
!                T3(thiscell) = T3(thiscell) +&
!                w2(thiscell, x) * c_E(3,1)*(dVxdx)+&
!                w2(thiscell, y) * c_E(3,2)*(dVydy)+&
!                w2(thiscell, z) * c_E(3,3)*(dVzdz) - e_piezo(3,3)*dEz(thiscell)*dt
!				
!                T4(thiscell) = T4(thiscell) +&
!                w2(thiscell, y) * c_E(4,4)*dVzdy+&
!                w2(thiscell, z) * c_E(4,4)*dVydz - &
!                dEy(thiscell)*e_piezo(2,4)*dt
!
!                T5(thiscell) = T5(thiscell) +&
!                w2(thiscell, x) * c_E(5,5)*dVzdx+&
!                w2(thiscell, z) * c_E(5,5)*dVxdz - &
!                dEx(thiscell)*e_piezo(1,5)*dt
!
!                T6(thiscell) = T6(thiscell) +&
!                w2(thiscell, x) * c_E(6,6)*dVydx+&
!                w2(thiscell, y) * c_E(6,6)*dVxdy
                
                ELSE
                !Cubic
                T1_x(ix, iy, iz) = T1_x(ix, iy, iz) * w1(ix, iy, iz, x) +&
                w2(ix, iy, iz, x) * c_E(1,1)*(dVxdx)
                T1_y(ix, iy, iz) = T1_y(ix, iy, iz) * w1(ix, iy, iz, y) +&
                w2(ix, iy, iz, y) * c_E(1,2)*(dVydy)
                T1_z(ix, iy, iz) = T1_z(ix, iy, iz) * w1(ix, iy, iz, z) +&
                w2(ix, iy, iz, z) * c_E(1,3)*(dVzdz)
                T1(ix, iy, iz) = T1_x(ix, iy, iz) + T1_y(ix, iy, iz) + T1_z(ix, iy, iz)

                T2_x(ix, iy, iz) = T2_x(ix, iy, iz) * w1(ix, iy, iz, x) +&
                w2(ix, iy, iz, x) * c_E(2,1)*(dVxdx)
                T2_y(ix, iy, iz) = T2_y(ix, iy, iz) * w1(ix, iy, iz, y) +&
                w2(ix, iy, iz, y) * c_E(2,2)*(dVydy)
                T2_z(ix, iy, iz) = T2_z(ix, iy, iz) * w1(ix, iy, iz, z) +&
                w2(ix, iy, iz, z) * c_E(2,3)*(dVzdz)
                T2(ix, iy, iz) = T2_x(ix, iy, iz) + T2_y(ix, iy, iz) + T2_z(ix, iy, iz)

                T3_x(ix, iy, iz) = T3_x(ix, iy, iz) * w1(ix, iy, iz, x) +&
                w2(ix, iy, iz, x) * c_E(3,1)*(dVxdx)
                T3_y(ix, iy, iz) = T3_y(ix, iy, iz) * w1(ix, iy, iz, y) +&
                w2(ix, iy, iz, y) * c_E(3,2)*(dVydy)
                T3_z(ix, iy, iz) = T3_z(ix, iy, iz) * w1(ix, iy, iz, z) +&
                w2(ix, iy, iz, z) * c_E(3,3)*(dVzdz)
                T3(ix, iy, iz) = T3_x(ix, iy, iz) + T3_y(ix, iy, iz) + T3_z(ix, iy, iz)
				
                T4_y(ix, iy, iz) = T4_y(ix, iy, iz) * w1(ix, iy, iz, y) +&
                w2(ix, iy, iz, y) * c_E(4,4)*dVzdy
                T4_z(ix, iy, iz) = T4_z(ix, iy, iz) * w1(ix, iy, iz, z) +&
                w2(ix, iy, iz, z) * c_E(4,4)*dVydz
                T4(ix, iy, iz) = T4_y(ix, iy, iz) + T4_z(ix, iy, iz) - &
                dEx(ix, iy, iz)*e_piezo(1,4)*dt

                T5_x(ix, iy, iz) = T5_x(ix, iy, iz) * w1(ix, iy, iz, x) +&
                w2(ix, iy, iz, x) * c_E(5,5)*dVzdx
                T5_z(ix, iy, iz) = T5_z(ix, iy, iz) * w1(ix, iy, iz, z) +&
                w2(ix, iy, iz, z) * c_E(5,5)*dVxdz
                T5(ix, iy, iz) = T5_x(ix, iy, iz) + T5_z(ix, iy, iz) - &
                dEy(ix, iy, iz)*e_piezo(2,5)*dt

                T6_x(ix, iy, iz) = T6_x(ix, iy, iz) * w1(ix, iy, iz, x) +&
                w2(ix, iy, iz, x) * c_E(6,6)*dVydx
                T6_y(ix, iy, iz) = T6_y(ix, iy, iz) * w1(ix, iy, iz, y) +&
                w2(ix, iy, iz, y) * c_E(6,6)*dVxdy
                T6(ix, iy, iz) = T6_x(ix, iy, iz) + T6_y(ix, iy, iz) - &
                dEz(ix, iy, iz)*e_piezo(3,6)*dt
                
                !Orthorhombic 2mm
!                T1_x(thiscell) = T1_x(thiscell) * w1(thiscell, x) +&
!                w2(thiscell, x) * c_E(1,1)*(dVxdx)
!                T1_y(thiscell) = T1_y(thiscell) * w1(thiscell, y) +&
!                w2(thiscell, y) * c_E(1,2)*(dVydy)
!                T1_z(thiscell) = T1_z(thiscell) * w1(thiscell, z) +&
!                w2(thiscell, z) * c_E(1,3)*(dVzdz)
!                T1(thiscell) = T1_x(thiscell) + T1_y(thiscell) + T1_z(thiscell) - e_piezo(3,1)*dEz(thiscell)*dt
!
!                T2_x(thiscell) = T2_x(thiscell) * w1(thiscell, x) +&
!                w2(thiscell, x) * c_E(2,1)*(dVxdx)
!                T2_y(thiscell) = T2_y(thiscell) * w1(thiscell, y) +&
!                w2(thiscell, y) * c_E(2,2)*(dVydy)
!                T2_z(thiscell) = T2_z(thiscell) * w1(thiscell, z) +&
!                w2(thiscell, z) * c_E(2,3)*(dVzdz)
!                T2(thiscell) = T2_x(thiscell) + T2_y(thiscell) + T2_z(thiscell) - e_piezo(3,2)*dEz(thiscell)*dt
!
!                T3_x(thiscell) = T3_x(thiscell) * w1(thiscell, x) +&
!                w2(thiscell, x) * c_E(3,1)*(dVxdx)
!                T3_y(thiscell) = T3_y(thiscell) * w1(thiscell, y) +&
!                w2(thiscell, y) * c_E(3,2)*(dVydy)
!                T3_z(thiscell) = T3_z(thiscell) * w1(thiscell, z) +&
!                w2(thiscell, z) * c_E(3,3)*(dVzdz)
!                T3(thiscell) = T3_x(thiscell) + T3_y(thiscell) + T3_z(thiscell) - e_piezo(3,3)*dEz(thiscell)*dt
!				
!                T4_y(thiscell) = T4_y(thiscell) * w1(thiscell, y) +&
!                w2(thiscell, y) * c_E(4,4)*dVzdy
!                T4_z(thiscell) = T4_z(thiscell) * w1(thiscell, z) +&
!                w2(thiscell, z) * c_E(4,4)*dVydz
!                T4(thiscell) = T4_y(thiscell) + T4_z(thiscell) - &
!                dEy(thiscell)*e_piezo(2,4)*dt
!
!                T5_x(thiscell) = T5_x(thiscell) * w1(thiscell, x) +&
!                w2(thiscell, x) * c_E(5,5)*dVzdx
!                T5_z(thiscell) = T5_z(thiscell) * w1(thiscell, z) +&
!                w2(thiscell, z) * c_E(5,5)*dVxdz
!                T5(thiscell) = T5_x(thiscell) + T5_z(thiscell) - &
!                dEx(thiscell)*e_piezo(1,5)*dt
!
!                T6_x(thiscell) = T6_x(thiscell) * w1(thiscell, x) +&
!                w2(thiscell, x) * c_E(6,6)*dVydx
!                T6_y(thiscell) = T6_y(thiscell) * w1(thiscell, y) +&
!                w2(thiscell, y) * c_E(6,6)*dVxdy
!                T6(thiscell) = T6_x(thiscell) + T6_y(thiscell)
                END IF
                
                Disx(ix, iy, iz)=Disx(ix, iy, iz)+ dDx(ix, iy, iz)
                Disy(ix, iy, iz)=Disy(ix, iy, iz)+ dDy(ix, iy, iz)
                Disz(ix, iy, iz)=Disz(ix, iy, iz)+ dDz(ix, iy, iz)
                
                Ex(ix, iy, iz)=Ex(ix, iy, iz)+ dEx(ix, iy, iz)
                Ey(ix, iy, iz)=Ey(ix, iy, iz)+ dEy(ix, iy, iz)
                Ez(ix, iy, iz)=Ez(ix, iy, iz)+ dEz(ix, iy, iz)

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
            T5(ix, iy, iz - 1) = -T5(ix, iy, iz)
            T3(ix, iy, iz) = 0
            T4(ix, iy, iz - 1) = -T4(ix, iy, iz)
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
        sendbuff_RIGHT(0*Nz*Ny+iz*Ny+iy)= Vx(Nx-2,iy,iz)
        sendbuff_RIGHT(1*Nz*Ny+iz*Ny+iy)= Vy(Nx-2,iy,iz)
        sendbuff_RIGHT(2*Nz*Ny+iz*Ny+iy)= Vz(Nx-2,iy,iz)
        sendbuff_LEFT(0*Nz*Ny+iz*Ny+iy) = Vx(1_li,iy,iz)
        sendbuff_LEFT(1*Nz*Ny+iz*Ny+iy) = Vy(1_li,iy,iz)
        sendbuff_LEFT(2*Nz*Ny+iz*Ny+iy) = Vz(1_li,iy,iz)
    END DO
    END DO
        
    DO iz=0, Nz-1
    DO ix=0, Nx-1
        sendbuff_UP(0*Nz*Nx+iz*Nx+ix)= Vx(ix,Ny-2,iz)
        sendbuff_UP(1*Nz*Nx+iz*Nx+ix)= Vy(ix,Ny-2,iz)
        sendbuff_UP(2*Nz*Nx+iz*Nx+ix)= Vz(ix,Ny-2,iz)
        sendbuff_DOWN(0*Nz*Nx+iz*Nx+ix) = Vx(ix,1_li,iz)
        sendbuff_DOWN(1*Nz*Nx+iz*Nx+ix) = Vy(ix,1_li,iz)
        sendbuff_DOWN(2*Nz*Nx+iz*Nx+ix) = Vz(ix,1_li,iz)
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
        sendbuff_RIGHT(0*Nz*Ny+iz*Ny+iy)= T1(Nx-2,iy,iz)
        sendbuff_RIGHT(1*Nz*Ny+iz*Ny+iy)= T2(Nx-2,iy,iz)
        sendbuff_RIGHT(2*Nz*Ny+iz*Ny+iy)= T3(Nx-2,iy,iz)
        sendbuff_RIGHT(3*Nz*Ny+iz*Ny+iy)= T4(Nx-2,iy,iz)
        sendbuff_RIGHT(4*Nz*Ny+iz*Ny+iy)= T5(Nx-2,iy,iz)
        sendbuff_RIGHT(5*Nz*Ny+iz*Ny+iy)= T6(Nx-2,iy,iz)
        sendbuff_LEFT(0*Nz*Ny+iz*Ny+iy) = T1(1_li,iy,iz)
        sendbuff_LEFT(1*Nz*Ny+iz*Ny+iy) = T2(1_li,iy,iz)
        sendbuff_LEFT(2*Nz*Ny+iz*Ny+iy) = T3(1_li,iy,iz)
        sendbuff_LEFT(3*Nz*Ny+iz*Ny+iy) = T4(1_li,iy,iz)
        sendbuff_LEFT(4*Nz*Ny+iz*Ny+iy) = T5(1_li,iy,iz)
        sendbuff_LEFT(5*Nz*Ny+iz*Ny+iy) = T6(1_li,iy,iz)
    END DO
    END DO
        
    DO iz=0, Nz-1
    DO ix=0, Nx-1
        sendbuff_UP(0*Nz*Nx+iz*Nx+ix)= T1(ix,Ny-2,iz)
        sendbuff_UP(1*Nz*Nx+iz*Nx+ix)= T2(ix,Ny-2,iz)
        sendbuff_UP(2*Nz*Nx+iz*Nx+ix)= T3(ix,Ny-2,iz)
        sendbuff_UP(3*Nz*Nx+iz*Nx+ix)= T4(ix,Ny-2,iz)
        sendbuff_UP(4*Nz*Nx+iz*Nx+ix)= T5(ix,Ny-2,iz)
        sendbuff_UP(5*Nz*Nx+iz*Nx+ix)= T6(ix,Ny-2,iz)
        sendbuff_DOWN(0*Nz*Nx+iz*Nx+ix) = T1(ix,1_li,iz)
        sendbuff_DOWN(1*Nz*Nx+iz*Nx+ix) = T2(ix,1_li,iz)
        sendbuff_DOWN(2*Nz*Nx+iz*Nx+ix) = T3(ix,1_li,iz)
        sendbuff_DOWN(3*Nz*Nx+iz*Nx+ix) = T4(ix,1_li,iz)
        sendbuff_DOWN(4*Nz*Nx+iz*Nx+ix) = T5(ix,1_li,iz)
        sendbuff_DOWN(5*Nz*Nx+iz*Nx+ix) = T6(ix,1_li,iz)
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






