!     
! File:   constants_module.f90
! Author: ludwig
!
! Created on April 12, 2014, 1:20 PM
!

MODULE constants_module
    USE Type_Kinds
    ! No implicit typing
    IMPLICIT NONE
    ! Explicit visibility declaration
    PRIVATE

    REAL(Double), PARAMETER, PUBLIC :: PI = 3.1415926535897932_dp
    REAL(Double), PARAMETER, PUBLIC :: clight = 299792458.0_dp
    COMPLEX(Double), PARAMETER, PUBLIC :: j1 = (0.0_dp, 1.0_dp)

    INTEGER(Short), PARAMETER, PUBLIC :: name_len = 25

    INTEGER, PARAMETER, PUBLIC :: x = 1
    INTEGER, PARAMETER, PUBLIC :: y = 2
    INTEGER, PARAMETER, PUBLIC :: z = 3
    
    !INTEGER, PARAMETER, PUBLIC :: Nstep=5000
    INTEGER, PARAMETER, PUBLIC :: Nprobe = 9*9*9
    
    CHARACTER(LEN=name_len), PARAMETER, PUBLIC :: PLM_archive='configfiles/PLM'
    CHARACTER(LEN=name_len), PARAMETER, PUBLIC :: sim_config='configfiles/sim.config'
    CHARACTER(LEN=10), PARAMETER, PUBLIC :: Materials_dir='Materials/'
    
    
    REAL(Double), PARAMETER, PUBLIC :: PWIDTH=0.4e-9	


END MODULE constants_module
