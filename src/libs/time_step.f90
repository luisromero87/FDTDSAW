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

    INTEGER(Long) UROLL3

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

    INTEGER(Long) UROLL3

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

    INTEGER(Long) UROLL3

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

    INTEGER(Long) UROLL3

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

    INTEGER(Long) UROLL3

    INTEGER :: ix, iy, iz

    ix = Nx/2
    iy = Ny/2
    iz = 1

    Vz(UROLL3(ix, iy, iz)) = Vz(UROLL3(ix, iy, iz)) + &
    (step * dt - 3 * PWIDTH)/(3 * PWIDTH) * exp(-1.0 * ((step * dt - 3.0 * PWIDTH)/(PWIDTH))**2)

    !WRITE(*, *) (step * dt - 3 * PWIDTH)/(3 * PWIDTH) * exp(-1.0 * ((step * dt - 3.0 * PWIDTH)/(PWIDTH))**2)

END SUBROUTINE dot_source