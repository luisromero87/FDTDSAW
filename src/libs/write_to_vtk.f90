!     
! File:   write_to_vtk.f90
! Author: ludwig
!
! Created on April 12, 2014, 9:35 PM
!

SUBROUTINE open_vtk_file(outfile)
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE
    
    CHARACTER(LEN=name_len), INTENT(IN) :: outfile
    
    INTEGER :: ix,iy,iz,axis
    
    2000 format(a)
    3000 format('DIMENSIONS ',i4,i4,i4)
    4000 format(1p,i4)
    
    OPEN(UNIT=12, FILE=outfile, ACTION="write", STATUS="replace")
    
    WRITE(12,2000) "# vtk DataFile Version 2.0"
    WRITE(12,2000) "Volume example"
    WRITE(12,2000) "ASCII"
    WRITE(12,2000) "DATASET RECTILINEAR_GRID"
    WRITE(12,3000) Nx,Ny-2,Nz-1
    
    WRITE(12,*) "X_COORDINATES", Nx, "double"
    DO ix = 0, Nx-1
      WRITE(12,4000) ix
    END DO
    WRITE(12,*) "Y_COORDINATES", Ny-2, "double"
    DO iy = 1, Ny-2
      WRITE(12,4000) iy
    END DO
    WRITE(12,*) "Z_COORDINATES", Nz-1, "double"
    DO iz = 1, Nz-1
      WRITE(12,4000) iz
    END DO
    
    CLOSE(12)
    
END SUBROUTINE open_vtk_file

SUBROUTINE write_point_data (outfile, data_name)
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE
    
    CHARACTER(LEN=name_len), INTENT(IN) :: data_name
    CHARACTER(LEN=name_len), INTENT(IN) :: outfile
    
    INTEGER(Long) UROLL3
    
    INTEGER :: ix,iy,iz,axis
    
    5000 format(1pe10.3,1x,1pe10.3,1x,1pe10.3)
    
    OPEN(UNIT=12, FILE=outfile, ACTION="write", STATUS="old", position="append")
    
    WRITE(12,*) "POINT_DATA", Nx*(Ny-2)*(Nz-1)
    WRITE(12,*) "SCALARS "//data_name//" double 3"
    WRITE(12,*) "LOOKUP_TABLE default"
    DO iz = 1, Nz-1
      DO iy = 1, Ny-2
	DO ix = 0, Nx-1
	    write(12,5000) REAL(Vx(UROLL3(ix,iy,iz)),Single),REAL(Vy(UROLL3(ix,iy,iz)),Single),REAL(Vz(UROLL3(ix,iy,iz)),Single)
	END DO
      END DO 
    END DO
    write(*,*)'file closed'
    
    CLOSE(12)
    
END SUBROUTINE write_point_data