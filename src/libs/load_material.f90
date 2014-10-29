!     
! File:   load_material.f90
! Author: ludwig
!
! Created on April 12, 2014, 9:35 PM
!

SUBROUTINE load_material ()
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE
    
    INTEGER :: st
            
    OPEN(UNIT = 10, FILE =  Materials_dir//material, FORM='unformatted',IOSTAT = st)
    READ(10) rho_inv
    READ(10) c_E
    READ(10) s_E
    READ(10) beta_s
    READ(10) e_piezo
    
!~     write (*,*) beta_s(1),beta_s(5),beta_s(9)
!~     write (*,*) e_piezo(10),e_piezo(14),e_piezo(12+6)
    
    CLOSE(10)
    
END SUBROUTINE load_material

