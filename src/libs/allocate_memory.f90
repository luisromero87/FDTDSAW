!     
! File:   allocate_memory.f90
! Author: ludwig
!
! Created on April 12, 2014, 9:35 PM
!

SUBROUTINE allocate_memory()
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE
    
    ALLOCATE(c_E(1:36))
    ALLOCATE(beta_s(1:9))
    ALLOCATE(e_piezo(1:18))
    
    ALLOCATE(w1(0:NCeldas - 1, 1:3))
    ALLOCATE(w2(0:NCeldas - 1, 1:3))
    
    ALLOCATE(dx(0:Nx - 1))
    ALLOCATE(dy(0:Ny - 1))
    ALLOCATE(dz(0:Nz - 1))

    ALLOCATE(Vx(0:NCeldas - 1))
    ALLOCATE(Vx_x(0:NCeldas - 1))
    ALLOCATE(Vx_y(0:NCeldas - 1))
    ALLOCATE(Vx_z(0:NCeldas - 1))

    ALLOCATE(Vy(0:NCeldas - 1))
    ALLOCATE(Vy_x(0:NCeldas - 1))
    ALLOCATE(Vy_y(0:NCeldas - 1))
    ALLOCATE(Vy_z(0:NCeldas - 1))

    ALLOCATE(Vz(0:NCeldas - 1))
    ALLOCATE(Vz_x(0:NCeldas - 1))
    ALLOCATE(Vz_y(0:NCeldas - 1))
    ALLOCATE(Vz_z(0:NCeldas - 1))

    ALLOCATE(T1(0:NCeldas - 1))
    ALLOCATE(T1_x(0:NCeldas - 1))
    ALLOCATE(T1_y(0:NCeldas - 1))
    ALLOCATE(T1_z(0:NCeldas - 1))

    ALLOCATE(T2(0:NCeldas - 1))
    ALLOCATE(T2_x(0:NCeldas - 1))
    ALLOCATE(T2_y(0:NCeldas - 1))
    ALLOCATE(T2_z(0:NCeldas - 1))

    ALLOCATE(T3(0:NCeldas - 1))
    ALLOCATE(T3_x(0:NCeldas - 1))
    ALLOCATE(T3_y(0:NCeldas - 1))
    ALLOCATE(T3_z(0:NCeldas - 1))

    ALLOCATE(T4(0:NCeldas - 1))
    ALLOCATE(T4_y(0:NCeldas - 1))
    ALLOCATE(T4_z(0:NCeldas - 1))

    ALLOCATE(T5(0:NCeldas - 1))
    ALLOCATE(T5_x(0:NCeldas - 1))
    ALLOCATE(T5_z(0:NCeldas - 1))

    ALLOCATE(T6(0:NCeldas - 1))
    ALLOCATE(T6_x(0:NCeldas - 1))
    ALLOCATE(T6_y(0:NCeldas - 1))


    
    c_E=0.0_dp
    beta_s=0.0_dp
    e_piezo=0.0_dp
    
    w1 = 1.0_dp
    w2 = 0.0_dp
    
    dx = 9.3333333_dp*1e-8
    dy = 9.3333333_dp*1e-8
    dz = 9.3333333_dp*1e-8
    
    Vx = 0.0_dp
    Vx_x = 0.0_dp
    Vx_y = 0.0_dp
    Vx_z = 0.0_dp

    Vy = 0.0_dp
    Vy_x = 0.0_dp
    Vy_y = 0.0_dp
    Vy_z = 0.0_dp

    Vz = 0.0_dp
    Vz_x = 0.0_dp
    Vz_y = 0.0_dp
    Vz_z = 0.0_dp

    T1 = 0.0_dp
    T1_x = 0.0_dp
    T1_y = 0.0_dp
    T1_z = 0.0_dp

    T2 = 0.0_dp
    T2_x = 0.0_dp
    T2_y = 0.0_dp
    T2_z = 0.0_dp

    T3 = 0.0_dp
    T3_x = 0.0_dp
    T3_y = 0.0_dp
    T3_z = 0.0_dp

    T4 = 0.0_dp
    T4_y = 0.0_dp
    T4_z = 0.0_dp

    T5 = 0.0_dp
    T5_x = 0.0_dp
    T5_z = 0.0_dp

    T6 = 0.0_dp
    T6_x = 0.0_dp
    T6_y = 0.0_dp


END SUBROUTINE allocate_memory

