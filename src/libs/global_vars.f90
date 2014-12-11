!     
! File:   global_vars.f90
! Author: ludwig
!
! Created on April 13, 2014, 2:11 PM
!

MODULE Global_Vars

#ifdef MPI2
    USE MPI
#endif
USE Type_Kinds
USE Constants_Module
IMPLICIT NONE


!MPI VARS
#ifdef MPI2
INTEGER, DIMENSION(MPI_STATUS_SIZE), PUBLIC :: mpistatus
INTEGER :: stats(MPI_STATUS_SIZE,8), reqs(8)
#endif
INTEGER, PUBLIC :: ierr, ntasks, me
INTEGER, PUBLIC :: nUP, nDOWN, nLEFT, nRIGHT
INTEGER(Long), PUBLIC :: vbuffsizex
INTEGER(Long), PUBLIC :: vbuffsizey
REAL(Double), DIMENSION (:), POINTER, PUBLIC :: mpibufferx
REAL(Double), DIMENSION (:), POINTER, PUBLIC :: mpibuffery

REAL(Double), DIMENSION (:), POINTER, PUBLIC :: recvbuff_UP
REAL(Double), DIMENSION (:), POINTER, PUBLIC :: recvbuff_DOWN
REAL(Double), DIMENSION (:), POINTER, PUBLIC :: recvbuff_RIGHT
REAL(Double), DIMENSION (:), POINTER, PUBLIC :: recvbuff_LEFT

REAL(Double), DIMENSION (:), POINTER, PUBLIC :: sendbuff_UP
REAL(Double), DIMENSION (:), POINTER, PUBLIC :: sendbuff_DOWN
REAL(Double), DIMENSION (:), POINTER, PUBLIC :: sendbuff_RIGHT
REAL(Double), DIMENSION (:), POINTER, PUBLIC :: sendbuff_LEFT


CHARACTER(LEN=name_len), PUBLIC :: Debug='False'

!INPUT PARAMETERS
CHARACTER(LEN=name_len), PUBLIC :: input_param_file='configfiles/simconfig'
CHARACTER(LEN=name_len), PUBLIC :: output_dir='outputdata/'
CHARACTER(LEN = name_len), PUBLIC :: material = 'bi12geo20' !Default
CHARACTER(LEN = name_len) :: D0_file = 'False' 

INTEGER(Long), PUBLIC :: NGx
INTEGER(Long), PUBLIC :: NGy
INTEGER(Long), PUBLIC :: NGz
INTEGER(Long), PUBLIC :: Nx
INTEGER(Long), PUBLIC :: Ny
INTEGER(Long), PUBLIC :: Nz
INTEGER(Long), PUBLIC :: NGCeldas
INTEGER(Long), PUBLIC :: NCeldas
INTEGER(Short), PUBLIC :: Nprocsx=1
INTEGER(Short), PUBLIC :: Nprocsy=1
INTEGER(Short), PUBLIC :: procsx=1
INTEGER(Short), PUBLIC :: procsy=1
REAL(Double), PUBLIC :: deltax
REAL(Double), PUBLIC :: deltay
REAL(Double), PUBLIC :: deltaz
INTEGER(Long), PUBLIC :: Nstep=5000
INTEGER(Long), PUBLIC :: data_fstep=100
REAL(Double), PUBLIC :: dt = 2.0_dp * 1.0e-12


!
INTEGER(Long), PUBLIC :: STEP
REAL(Double), PUBLIC :: offsetx
REAL(Double), PUBLIC :: offsety
REAL(Double), DIMENSION(:), POINTER, PUBLIC :: dx
REAL(Double), DIMENSION(:), POINTER, PUBLIC :: dy
REAL(Double), DIMENSION(:), POINTER, PUBLIC :: dz

!Material properties
REAL(Double), PUBLIC :: rho, rho_inv = 0.0_dp
REAL(Double), DIMENSION(:,:), POINTER, PUBLIC :: c_E, s_E
REAL(Double), DIMENSION(:,:), POINTER, PUBLIC :: beta_s
REAL(Double), DIMENSION(:,:), POINTER, PUBLIC :: e_piezo

!Mechanical variables
REAL(Double), DIMENSION(:,:,:), POINTER, PUBLIC :: Vx
REAL(Double), DIMENSION(:,:,:), POINTER, PUBLIC :: Vx_x
REAL(Double), DIMENSION(:,:,:), POINTER, PUBLIC :: Vx_y
REAL(Double), DIMENSION(:,:,:), POINTER, PUBLIC :: Vx_z

REAL(Double), DIMENSION(:,:,:), POINTER, PUBLIC :: Vy
REAL(Double), DIMENSION(:,:,:), POINTER, PUBLIC :: Vy_x
REAL(Double), DIMENSION(:,:,:), POINTER, PUBLIC :: Vy_y
REAL(Double), DIMENSION(:,:,:), POINTER, PUBLIC :: Vy_z

REAL(Double), DIMENSION(:,:,:), POINTER, PUBLIC :: Vz
REAL(Double), DIMENSION(:,:,:), POINTER, PUBLIC :: Vz_x
REAL(Double), DIMENSION(:,:,:), POINTER, PUBLIC :: Vz_y
REAL(Double), DIMENSION(:,:,:), POINTER, PUBLIC :: Vz_z

REAL(Double), DIMENSION(:,:,:), POINTER :: T1
REAL(Double), DIMENSION(:,:,:), POINTER :: T1_x
REAL(Double), DIMENSION(:,:,:), POINTER :: T1_y
REAL(Double), DIMENSION(:,:,:), POINTER :: T1_z

REAL(Double), DIMENSION(:,:,:), POINTER :: T2
REAL(Double), DIMENSION(:,:,:), POINTER :: T2_x
REAL(Double), DIMENSION(:,:,:), POINTER :: T2_y
REAL(Double), DIMENSION(:,:,:), POINTER :: T2_z

REAL(Double), DIMENSION(:,:,:), POINTER :: T3
REAL(Double), DIMENSION(:,:,:), POINTER :: T3_x
REAL(Double), DIMENSION(:,:,:), POINTER :: T3_y
REAL(Double), DIMENSION(:,:,:), POINTER :: T3_z

REAL(Double), DIMENSION(:,:,:), POINTER :: T4
REAL(Double), DIMENSION(:,:,:), POINTER :: T4_y
REAL(Double), DIMENSION(:,:,:), POINTER :: T4_z

REAL(Double), DIMENSION(:,:,:), POINTER :: T5
REAL(Double), DIMENSION(:,:,:), POINTER :: T5_x
REAL(Double), DIMENSION(:,:,:), POINTER :: T5_z

REAL(Double), DIMENSION(:,:,:), POINTER :: T6
REAL(Double), DIMENSION(:,:,:), POINTER :: T6_x
REAL(Double), DIMENSION(:,:,:), POINTER :: T6_y

!Electrical variable
REAL(Double), DIMENSION(:,:,:), POINTER :: Ex
REAL(Double), DIMENSION(:,:,:), POINTER :: Ey
REAL(Double), DIMENSION(:,:,:), POINTER :: Ez

REAL(Double), DIMENSION(:,:,:), POINTER :: dEx
REAL(Double), DIMENSION(:,:,:), POINTER :: dEy
REAL(Double), DIMENSION(:,:,:), POINTER :: dEz

REAL(Double), DIMENSION(:,:,:), POINTER :: dDx
REAL(Double), DIMENSION(:,:,:), POINTER :: dDy
REAL(Double), DIMENSION(:,:,:), POINTER :: dDz

REAL(Double), DIMENSION(:,:,:), ALLOCATABLE :: D0x
REAL(Double), DIMENSION(:,:,:), ALLOCATABLE :: D0y
REAL(Double), DIMENSION(:,:,:), ALLOCATABLE :: D0z

REAL(Double), DIMENSION(:,:,:), POINTER :: Disx
REAL(Double), DIMENSION(:,:,:), POINTER :: Disy
REAL(Double), DIMENSION(:,:,:), POINTER :: Disz

REAL(Double) :: phase=0.0_dp

!PML vars

INTEGER :: PMLwidth
        
REAL(Double) :: m, smax

INTEGER :: xstart,xend,ystart,yend,zstart,zend

REAL(Double), DIMENSION(:,:,:,:), POINTER, PUBLIC :: w1
REAL(Double), DIMENSION(:,:,:,:), POINTER, PUBLIC :: w2

REAL(Double), PUBLIC :: U_k
REAL(Double), PUBLIC :: U_k_total

REAL(Double), PUBLIC :: S1,S2,S3,S4,S5,S6
REAL(Double), PUBLIC :: U_s
REAL(Double), PUBLIC :: U_s_total

REAL(Double), PUBLIC :: U_e
REAL(Double), PUBLIC :: U_e_total

ENDMODULE Global_Vars

