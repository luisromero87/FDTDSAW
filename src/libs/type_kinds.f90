!     
! File:   type_kinds.f90
! Author: ludwig
!
! Created on April 12, 2014, 1:20 PM
!

MODULE Type_Kinds
    ! No implicit typing
    IMPLICIT NONE
    ! Explicit visibility declaration
    PRIVATE
    PUBLIC :: Byte, Short, Long
    PUBLIC :: Single, Double
    PUBLIC :: li, dp
    ! Integer kinds
    INTEGER, PARAMETER :: Byte = SELECTED_INT_KIND(1) ! Byte
    INTEGER, PARAMETER :: Short = SELECTED_INT_KIND(4) ! Short
    INTEGER, PARAMETER :: Long = SELECTED_INT_KIND(8) ! Long
    ! Floating point kinds
    INTEGER, PARAMETER :: Single = SELECTED_REAL_KIND(6) ! Single
    INTEGER, PARAMETER :: Double = SELECTED_REAL_KIND(15) ! Double
    ! Generic kinds
    INTEGER, PARAMETER :: li = Long ! Generic integer kind
    INTEGER, PARAMETER :: dp = Double ! Generic real kind
    
END MODULE Type_Kinds
