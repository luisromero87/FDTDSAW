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
    4000 format(1p,e12.6)
    
    OPEN(UNIT=12, FILE='outputdata/'//outfile, ACTION="write", STATUS="replace")
    
    WRITE(12,2000) "# vtk DataFile Version 2.0"
    WRITE(12,2000) "Volume data"
    WRITE(12,2000) "ASCII"
    WRITE(12,2000) "DATASET RECTILINEAR_GRID"
    WRITE(12,3000) Nx-2,Ny-2,Nz-1
    
    WRITE(12,*) "X_COORDINATES", Nx-2, "double"
    DO ix = 0, Nx-3
      WRITE(12,4000) ix*deltax+offsetx
    END DO
    WRITE(12,*) "Y_COORDINATES", Ny-2, "double"
    DO iy = 0, Ny-3
      WRITE(12,4000) iy*deltay+offsety
    END DO
    WRITE(12,*) "Z_COORDINATES", Nz-1, "double"
    DO iz = 1, Nz-1
      WRITE(12,4000) iz*deltaz
    END DO
    
    CLOSE(12)
    
END SUBROUTINE open_vtk_file

SUBROUTINE write_volume_v (outfile, data_name)
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE
    
    CHARACTER(LEN=name_len), INTENT(IN) :: data_name
    CHARACTER(LEN=name_len), INTENT(IN) :: outfile
    
    INTEGER(Long) UROLL3
    
    INTEGER :: ix,iy,iz,axis
    
    5000 format(1pe10.3,1x,1pe10.3,1x,1pe10.3)
    
    OPEN(UNIT=12, FILE='outputdata/'//outfile, ACTION="write", STATUS="old", position="append")
    
    WRITE(12,*) "POINT_DATA", (Nx-2)*(Ny-2)*(Nz-1)
    WRITE(12,*) "SCALARS "//"v"//" double 3"
    WRITE(12,*) "LOOKUP_TABLE default"
    DO iz = 1, Nz-1
      DO iy = 1, Ny-2
	DO ix = 1, Nx-2
	    write(12,5000) REAL(Vx(ix,iy,iz),Single),REAL(Vy(ix,iy,iz),Single),REAL(Vz(ix,iy,iz),Single)
	END DO
      END DO 
    END DO
    
    WRITE(12,*) "SCALARS "//"E"//" double 3"
    WRITE(12,*) "LOOKUP_TABLE default"
    DO iz = 1, Nz-1
      DO iy = 1, Ny-2
	DO ix = 1, Nx-2
	    write(12,5000) REAL(Ex(ix,iy,iz),Single),REAL(Ey(ix,iy,iz),Single),REAL(Ez(ix,iy,iz),Single)
	END DO
      END DO 
    END DO
!~     write(*,*)'file closed'
    
    CLOSE(12)
    
END SUBROUTINE write_volume_v

SUBROUTINE write_free_surface (outfile, data_name)
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE
    
    CHARACTER(LEN=name_len), INTENT(IN) :: outfile
    CHARACTER(LEN=name_len), INTENT(IN) :: data_name
    
    INTEGER(Long) UROLL3
    
    INTEGER :: ix,iy,iz,axis
    
    2000 format(a)
    3000 format('DIMENSIONS ',i4,i4,i4)
    4000 format(1p,e12.6)
    5000 format(1pe10.3,1x,1pe10.3,1x,1pe10.3)
    
    
    OPEN(UNIT=12, FILE='outputdata/'//outfile, ACTION="write", STATUS="replace")
    
    WRITE(12,2000) "# vtk DataFile Version 2.0"
    WRITE(12,2000) "Free Surface"
    WRITE(12,2000) "ASCII"
    WRITE(12,2000) "DATASET RECTILINEAR_GRID"
    WRITE(12,3000) Nx-2,Ny-2,1
    
    WRITE(12,*) "X_COORDINATES", Nx-2, "double"
    DO ix = 0, Nx-3
      WRITE(12,4000) ix*deltax+offsetx
    END DO
    WRITE(12,*) "Y_COORDINATES", Ny-2, "double"
    DO iy = 0, Ny-3
      WRITE(12,4000) iy*deltay+offsety
    END DO
    WRITE(12,*) "Z_COORDINATES", 1, "double"
    DO iz = Nz/2, Nz/2
      WRITE(12,4000) iz*deltaz
    END DO
    
    CLOSE(12)
    
    
    OPEN(UNIT=12, FILE='outputdata/'//outfile, ACTION="write", STATUS="old", position="append")
    
    WRITE(12,*) "POINT_DATA", (Nx-2)*(Ny-2)*(1)
    WRITE(12,*) "SCALARS "//data_name//" double 3"
    WRITE(12,*) "LOOKUP_TABLE default"
    iz = Nz/2
    DO iy = 1, Ny-2
	DO ix = 1, Nx-2
	    write(12,5000) REAL(Vx(ix,iy,iz),Single),REAL(Vy(ix,iy,iz),Single),REAL(Vz(ix,iy,iz),Single)
	END DO
    END DO 
    
!~     write(*,*)'file closed'
    
    CLOSE(12)
    
END SUBROUTINE write_free_surface

SUBROUTINE write_volume_w1 ()
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE
    
    INTEGER(Long) UROLL3
    
    CHARACTER(LEN=name_len) :: data_name ='w1'
    CHARACTER(LEN=name_len) :: outfile
    
    INTEGER :: ix,iy,iz,axis
    
    WRITE(outfile, '(A,i3.3,A)') 'w1_',me,'.vtk' 
    
    2000 format(a)
    3000 format('DIMENSIONS ',i4,i4,i4)
    4000 format(1p,e12.6)
    
    OPEN(UNIT=12, FILE='outputdata/'//outfile, ACTION="write", STATUS="replace")
    
    WRITE(12,2000) "# vtk DataFile Version 2.0"
    WRITE(12,2000) "Volume data"
    WRITE(12,2000) "ASCII"
    WRITE(12,2000) "DATASET RECTILINEAR_GRID"
    WRITE(12,3000) Nx-2,Ny-2,Nz-1
    
    WRITE(12,*) "X_COORDINATES", Nx-2, "double"
    DO ix = 0, Nx-3
      WRITE(12,4000) ix*deltax+offsetx
    END DO
    WRITE(12,*) "Y_COORDINATES", Ny-2, "double"
    DO iy = 0, Ny-3
      WRITE(12,4000) iy*deltay+offsety
    END DO
    WRITE(12,*) "Z_COORDINATES", Nz-1, "double"
    DO iz = 1, Nz-1
      WRITE(12,4000) iz*deltaz
    END DO
    
    5000 format(1pe10.3,1x,1pe10.3,1x,1pe10.3)
        
    WRITE(12,*) "POINT_DATA", (Nx-2)*(Ny-2)*(Nz-1)
    WRITE(12,*) "SCALARS "//data_name//" double 3"
    WRITE(12,*) "LOOKUP_TABLE default"
    DO iz = 1, Nz-1
    DO iy = 1, Ny-2
	DO ix = 1, Nx-2
	    write(12,5000) REAL(w1(ix,iy,iz,x),Single),REAL(w1(ix,iy,iz,y),Single),REAL(w1(ix,iy,iz,z),Single)
	END DO
    END DO 
    END DO
    
    CLOSE(12)
    
END SUBROUTINE write_volume_w1

SUBROUTINE xzplane()
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE
    
    INTEGER(Long) UROLL3
    
    CHARACTER(LEN=name_len) :: data_name ='v'
    CHARACTER(LEN=name_len) :: outfile
    
    INTEGER :: ix,iy,iz,axis
    
    
	IF (procsy.EQ.(Nprocsy-1)/2) THEN
	    iy = (Ny/2-2)*MOD(Nprocsy-1,2)+Ny/2
    
        WRITE(outfile, '(A,i4.4,A,i4.4,A)') 'v_yzplane_',me,'_',step/10,'.vtk' 
        
        2000 format(a)
        3000 format('DIMENSIONS ',i4,i4,i4)
        4000 format(1p,e12.6)
        
        OPEN(UNIT=12, FILE='outputdata/'//outfile, ACTION="write", STATUS="replace")
        
        WRITE(12,2000) "# vtk DataFile Version 2.0"
        WRITE(12,2000) "Volume data"
        WRITE(12,2000) "ASCII"
        WRITE(12,2000) "DATASET RECTILINEAR_GRID"
        WRITE(12,3000) Nx-2,1,Nz-1
        
        WRITE(12,*) "X_COORDINATES", Nx-2, "double"
        DO ix = 0, Nx-3
          WRITE(12,4000) ix*deltax+offsetx
        END DO
        WRITE(12,*) "Y_COORDINATES", 1, "double"
!        DO iy = 0, Ny-3
          WRITE(12,4000) iy*deltay+offsety
!        END DO
        WRITE(12,*) "Z_COORDINATES", Nz-1, "double"
        DO iz = 1, Nz-1
          WRITE(12,4000) iz*deltaz
        END DO
        
        5000 format(1pe10.3,1x,1pe10.3,1x,1pe10.3)
            
        WRITE(12,*) "POINT_DATA", (Nx-2)*(1)*(Nz-1)
        WRITE(12,*) "SCALARS "//"v"//" double 3"
        WRITE(12,*) "LOOKUP_TABLE default"
                
        DO iz = 1, Nz-1
        DO ix = 1, Nx-2
            write(12,5000) REAL(Vx(ix,iy,iz),Single),REAL(Vy(ix,iy,iz),Single),REAL(Vz(ix,iy,iz),Single)
        END DO 
        END DO
        WRITE(12,*) "SCALARS "//"E"//" double 3"
        WRITE(12,*) "LOOKUP_TABLE default"
                
        DO iz = 1, Nz-1
        DO ix = 1, Nx-2
            write(12,5000) REAL(Ex(ix,iy,iz),Single),REAL(Ey(ix,iy,iz),Single),REAL(Ez(ix,iy,iz),Single)
        END DO 
        END DO
        
        CLOSE(12)
    END IF
    
END SUBROUTINE xzplane

SUBROUTINE write_volume_D0 ()
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE
    
    INTEGER(Long) UROLL3
    
    CHARACTER(LEN=name_len) :: data_name ='D0'
    CHARACTER(LEN=name_len) :: outfile
    
    INTEGER :: ix,iy,iz,axis
    
    WRITE(outfile, '(A,i3.3,A)') 'D0_',me,'.vtk' 
    
    2000 format(a)
    3000 format('DIMENSIONS ',i4,i4,i4)
    4000 format(1p,e12.6)
    
    OPEN(UNIT=12, FILE='outputdata/IDT/'//outfile, ACTION="write", STATUS="replace")
    
    WRITE(12,2000) "# vtk DataFile Version 2.0"
    WRITE(12,2000) "Volume data"
    WRITE(12,2000) "ASCII"
    WRITE(12,2000) "DATASET RECTILINEAR_GRID"
    WRITE(12,3000) Nx-2,Ny-2,Nz-1
    
    WRITE(12,*) "X_COORDINATES", Nx-2, "double"
    DO ix = 0, Nx-3
      WRITE(12,4000) ix*deltax+offsetx
    END DO
    WRITE(12,*) "Y_COORDINATES", Ny-2, "double"
    DO iy = 0, Ny-3
      WRITE(12,4000) iy*deltay+offsety
    END DO
    WRITE(12,*) "Z_COORDINATES", Nz-1, "double"
    DO iz = 1, Nz-1
      WRITE(12,4000) iz*deltaz
    END DO
    
    5000 format(1pe10.3,1x,1pe10.3,1x,1pe10.3)
        
    WRITE(12,*) "POINT_DATA", (Nx-2)*(Ny-2)*(Nz-1)
    WRITE(12,*) "SCALARS "//data_name//" double 3"
    WRITE(12,*) "LOOKUP_TABLE default"
    DO iz = 1, Nz-1
    DO iy = 1, Ny-2
	DO ix = 1, Nx-2
	    write(12,5000) REAL(D0x(ix,iy,iz),Single),REAL(D0y(ix,iy,iz),Single),REAL(D0z(ix,iy,iz),Single)
	END DO
    END DO 
    END DO
    
    CLOSE(12)
    
END SUBROUTINE write_volume_D0
