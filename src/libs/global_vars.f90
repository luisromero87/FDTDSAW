!     
! File:   global_vars.f90
! Author: ludwig
!
! Created on April 13, 2014, 2:11 PM
!

MODULE global_vars

    USE MPI
    USE Type_Kinds
    USE Constants_Module
    IMPLICIT NONE
    
    
!~     REAL(Double), DIMENSION(:), POINTER :: mpibufferx, mpibuffery
    INTEGER, DIMENSION(MPI_STATUS_SIZE), PUBLIC :: mpistatus
    INTEGER, PUBLIC :: ierr, ntasks, me, req1, req2
    
    CHARACTER(LEN = name_len), PUBLIC :: material = 'bi12geo20' !Default

    INTEGER(Long), PUBLIC :: Nx
    INTEGER(Long), PUBLIC :: Ny
    INTEGER(Long), PUBLIC :: Nz
    INTEGER(Long), PUBLIC :: NCeldas
    INTEGER(Short), PUBLIC :: Nprocsx=1
    INTEGER(Short), PUBLIC :: Nprocsy=1
    INTEGER(Short), PUBLIC :: procsx=1
    INTEGER(Short), PUBLIC :: procsy=1
    INTEGER(Short), PUBLIC :: nUP, nDOWN, nLEFT, nRIGHT

    INTEGER(Long), PUBLIC :: STEP
    INTEGER(Long), PUBLIC :: Nstep=5000
    REAL(Double), PUBLIC :: dt = 2.0_dp * 1.0e-12
    REAL(Double), PUBLIC :: deltax
    REAL(Double), PUBLIC :: deltay
    REAL(Double), PUBLIC :: deltaz
    REAL(Double), PUBLIC :: offsetx
    REAL(Double), PUBLIC :: offsety
    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: dx
    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: dy
    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: dz

    !Material properties
    REAL(Double), PUBLIC :: rho_inv = 0.0_dp
    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: c_E
    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: beta_s
    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: e_piezo

    !Mechanical variables
    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: Vx
    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: Vx_x
    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: Vx_y
    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: Vx_z

    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: Vy
    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: Vy_x
    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: Vy_y
    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: Vy_z

    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: Vz
    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: Vz_x
    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: Vz_y
    REAL(Double), DIMENSION(:), POINTER, PUBLIC :: Vz_z

    REAL(Double), DIMENSION(:), POINTER :: T1
    REAL(Double), DIMENSION(:), POINTER :: T1_x
    REAL(Double), DIMENSION(:), POINTER :: T1_y
    REAL(Double), DIMENSION(:), POINTER :: T1_z

    REAL(Double), DIMENSION(:), POINTER :: T2
    REAL(Double), DIMENSION(:), POINTER :: T2_x
    REAL(Double), DIMENSION(:), POINTER :: T2_y
    REAL(Double), DIMENSION(:), POINTER :: T2_z

    REAL(Double), DIMENSION(:), POINTER :: T3
    REAL(Double), DIMENSION(:), POINTER :: T3_x
    REAL(Double), DIMENSION(:), POINTER :: T3_y
    REAL(Double), DIMENSION(:), POINTER :: T3_z

    REAL(Double), DIMENSION(:), POINTER :: T4
    REAL(Double), DIMENSION(:), POINTER :: T4_y
    REAL(Double), DIMENSION(:), POINTER :: T4_z

    REAL(Double), DIMENSION(:), POINTER :: T5
    REAL(Double), DIMENSION(:), POINTER :: T5_x
    REAL(Double), DIMENSION(:), POINTER :: T5_z

    REAL(Double), DIMENSION(:), POINTER :: T6
    REAL(Double), DIMENSION(:), POINTER :: T6_x
    REAL(Double), DIMENSION(:), POINTER :: T6_y
    
    !Electrical variable
    REAL(Double), DIMENSION(:), POINTER :: dEx
    REAL(Double), DIMENSION(:), POINTER :: dEy
    REAL(Double), DIMENSION(:), POINTER :: dEz

    REAL(Double), DIMENSION(:), POINTER :: dDx
    REAL(Double), DIMENSION(:), POINTER :: dDy
    REAL(Double), DIMENSION(:), POINTER :: dDz

    REAL(Double), DIMENSION(:), POINTER :: D0x
    REAL(Double), DIMENSION(:), POINTER :: D0y
    REAL(Double), DIMENSION(:), POINTER :: D0z
    
    REAL(Double) :: phase=0.0_dp

    !PML weights
    REAL(Double), DIMENSION(:,:), POINTER, PUBLIC :: w1
    REAL(Double), DIMENSION(:,:), POINTER, PUBLIC :: w2

END MODULE global_vars

