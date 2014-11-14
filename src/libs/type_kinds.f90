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
PUBLIC :: FR8P, FR4P, FR_P
PUBLIC :: DR8P, DR4P, DR_P
PUBLIC ::   FI8P,   FI4P,   FI2P,   FI1P,   FI_P
PUBLIC :: FI8PZP, FI4PZP, FI2PZP, FI1PZP, FI_PZP
PUBLIC ::   DI8P,   DI4P,   DI2P,   DI1P,   DI_P

! Integer kinds
INTEGER, PARAMETER :: Byte = SELECTED_INT_KIND(1) ! Byte
INTEGER, PARAMETER :: Short = SELECTED_INT_KIND(4) ! Short
INTEGER, PARAMETER :: Long = SELECTED_INT_KIND(8) ! Long

! Floating point kinds
INTEGER, PARAMETER :: Single = SELECTED_REAL_KIND(6,37) !< 6   digits, range \f$[10^{-37}  , 10^{+37}   - 1]\f$; 32 bits.
INTEGER, PARAMETER :: Double = SELECTED_REAL_KIND(15,307) !< 15  digits, range \f$[10^{-307} , 10^{+307}  - 1]\f$; 64 bits.

! Generic kinds
INTEGER, PARAMETER :: li = Long ! Generic integer kind
INTEGER, PARAMETER :: dp = Double ! Generic real kind

! Format parameters useful for writing in a well-ascii-format numeric variables.
! Real output formats:
character(10), parameter:: FR8P  = '(E23.15E3)' !< Output format for kind=R8P variable.
character(9),  parameter:: FR4P  = '(E13.6E2)'  !< Output format for kind=R4P variable.
character(10), parameter:: FR_P  = FR8P         !< Output format for kind=R_P variable.
! Real number of digits of output formats:
integer, parameter:: DR8P  = 23   !< Number of digits of output format FR8P.
integer, parameter:: DR4P  = 13   !< Number of digits of output format FR4P.
integer, parameter:: DR_P  = DR8P !< Number of digits of output format FR_P.
! Integer output formats:
character(5), parameter:: FI8P   = '(I20)'    !< Output format                     for kind=I8P variable.
character(8), parameter:: FI8PZP = '(I20.19)' !< Output format with zero prefixing for kind=I8P variable.
character(5), parameter:: FI4P   = '(I11)'    !< Output format                     for kind=I4P variable.
character(8), parameter:: FI4PZP = '(I11.10)' !< Output format with zero prefixing for kind=I4P variable.
character(4), parameter:: FI2P   = '(I6)'     !< Output format                     for kind=I2P variable.
character(6), parameter:: FI2PZP = '(I6.5)'   !< Output format with zero prefixing for kind=I2P variable.
character(4), parameter:: FI1P   = '(I4)'     !< Output format                     for kind=I1P variable.
character(6), parameter:: FI1PZP = '(I4.3)'   !< Output format with zero prefixing for kind=I1P variable.
character(5), parameter:: FI_P   = FI4P       !< Output format                     for kind=I_P variable.
character(8), parameter:: FI_PZP = FI4PZP     !< Output format with zero prefixing for kind=I_P variable.
! Integer number of digits of output formats:
integer, parameter:: DI8P = 20   !< Number of digits of output format I8P.
integer, parameter:: DI4P = 11   !< Number of digits of output format I4P.
integer, parameter:: DI2P = 6    !< Number of digits of output format I2P.
integer, parameter:: DI1P = 4    !< Number of digits of output format I1P.
integer, parameter:: DI_P = DI4P !< Number of digits of output format I_P.
    
END MODULE Type_Kinds
