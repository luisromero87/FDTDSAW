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

    ALLOCATE(dEx(0:NCeldas - 1))
    ALLOCATE(dEy(0:NCeldas - 1))
    ALLOCATE(dEz(0:NCeldas - 1))

    ALLOCATE(dDx(0:NCeldas - 1))
    ALLOCATE(dDy(0:NCeldas - 1))
    ALLOCATE(dDz(0:NCeldas - 1))

    ALLOCATE(D0x(0:NCeldas - 1))
    ALLOCATE(D0y(0:NCeldas - 1))
    ALLOCATE(D0z(0:NCeldas - 1))
    
    ALLOCATE(w1(0:NCeldas - 1, 1:3))
    ALLOCATE(w2(0:NCeldas - 1, 1:3))
    
    ALLOCATE(mpibufferx(0:6*Nz*Ny-1))
    ALLOCATE(mpibuffery(0:6*Nz*Nx-1))
    
    ALLOCATE(recvbuff_UP(0:6*Nz*Nx-1))
    ALLOCATE(recvbuff_DOWN(0:6*Nz*Nx-1))
    ALLOCATE(recvbuff_RIGHT(0:6*Nz*Ny-1))
    ALLOCATE(recvbuff_LEFT(0:6*Nz*Ny-1))
    
    ALLOCATE(sendbuff_UP(0:6*Nz*Nx-1))
    ALLOCATE(sendbuff_DOWN(0:6*Nz*Nx-1))
    ALLOCATE(sendbuff_RIGHT(0:6*Nz*Ny-1))
    ALLOCATE(sendbuff_LEFT(0:6*Nz*Ny-1))


    
    c_E=0.0_dp
    beta_s=0.0_dp
    e_piezo=0.0_dp
    
    w1 = 1.0_dp
    w2 = 0.0_dp
    
    dx = deltax
    dy = deltay
    dz = deltaz
    
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

    dEx = 0.0_dp
    dEy = 0.0_dp
    dEz = 0.0_dp

    dDx = 0.0_dp
    dDy = 0.0_dp
    dDz = 0.0_dp

    D0x = 0.0_dp
    D0y = 0.0_dp
    D0z = 0.0_dp


END SUBROUTINE allocate_memory

SUBROUTINE deallocate_memory()
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE
    
!~     ALLOCATE(mpibufferx(0:3*Ny*Nz-1))
!~     ALLOCATE(mpibuffery(0:Nx*Nz-1))
    
    DEALLOCATE(c_E)
    DEALLOCATE(beta_s)
    DEALLOCATE(e_piezo)
    
    DEALLOCATE(dx)
    DEALLOCATE(dy)
    DEALLOCATE(dz)

    DEALLOCATE(Vx)
    DEALLOCATE(Vx_x)
    DEALLOCATE(Vx_y)
    DEALLOCATE(Vx_z)

    DEALLOCATE(Vy)
    DEALLOCATE(Vy_x)
    DEALLOCATE(Vy_y)
    DEALLOCATE(Vy_z)

    DEALLOCATE(Vz)
    DEALLOCATE(Vz_x)
    DEALLOCATE(Vz_y)
    DEALLOCATE(Vz_z)

    DEALLOCATE(T1)
    DEALLOCATE(T1_x)
    DEALLOCATE(T1_y)
    DEALLOCATE(T1_z)

    DEALLOCATE(T2)
    DEALLOCATE(T2_x)
    DEALLOCATE(T2_y)
    DEALLOCATE(T2_z)

    DEALLOCATE(T3)
    DEALLOCATE(T3_x)
    DEALLOCATE(T3_y)
    DEALLOCATE(T3_z)

    DEALLOCATE(T4)
    DEALLOCATE(T4_y)
    DEALLOCATE(T4_z)

    DEALLOCATE(T5)
    DEALLOCATE(T5_x)
    DEALLOCATE(T5_z)

    DEALLOCATE(T6)
    DEALLOCATE(T6_x)
    DEALLOCATE(T6_y)

    DEALLOCATE(dEx)
    DEALLOCATE(dEy)
    DEALLOCATE(dEz)

    DEALLOCATE(dDx)
    DEALLOCATE(dDy)
    DEALLOCATE(dDz)

    DEALLOCATE(D0x)
    DEALLOCATE(D0y)
    DEALLOCATE(D0z)

    DEALLOCATE(w1)
    DEALLOCATE(w2)

END SUBROUTINE deallocate_memory
