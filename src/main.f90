PROGRAM acousticwaves

    !USE mpi
    USE Type_Kinds
    USE Constants_Module
    USE Global_Vars
    IMPLICIT NONE

    INTERFACE

        SUBROUTINE load_material(material, rho_inv, c_E, beta_s, e_piezo)
            USE Type_Kinds
            USE Constants_Module
            IMPLICIT NONE
            CHARACTER(LEN = name_len) :: material
            REAL(Double) :: rho_inv
            REAL(Double), DIMENSION(:), POINTER :: c_E
            REAL(Double), DIMENSION(:), POINTER :: beta_s
            REAL(Double), DIMENSION(:), POINTER :: e_piezo
        END SUBROUTINE load_material

        SUBROUTINE load_PML()
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
        END SUBROUTINE load_PML

        SUBROUTINE set_up_sim(material, NCeldas, Nx, Ny, Nz)
            USE Type_Kinds
            USE Constants_Module
            IMPLICIT NONE
            CHARACTER(LEN = name_len) :: material
            INTEGER(Long) :: NCeldas
            INTEGER(Long) :: Nx
            INTEGER(Long) :: Ny
            INTEGER(Long) :: Nz
        END SUBROUTINE set_up_sim

        SUBROUTINE allocate_memory()
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
        END SUBROUTINE allocate_memory

        SUBROUTINE open_vtk_file(outfile)
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE

            CHARACTER(LEN = name_len), INTENT(IN) :: outfile
        END SUBROUTINE

        SUBROUTINE write_point_data(outfile, data_name)
            USE Type_Kinds
            USE Constants_Module
            USE Global_Vars
            IMPLICIT NONE
            CHARACTER(LEN = name_len), INTENT(IN) :: outfile
            CHARACTER(LEN = name_len), INTENT(IN) :: data_name
        END SUBROUTINE
    END INTERFACE

    !Functions
    INTEGER(Long) :: UROLL3

    CHARACTER(LEN = name_len) :: material = 'bi12geo20' !Default
    CHARACTER(LEN = name_len) :: outfile = 'prueba.vtk' !Default
    CHARACTER(LEN = name_len) :: data_name = 'w1' !Default

    INTEGER :: st

    INTEGER :: axis, index, ix, iy, iz, cont

    REAL(Double), DIMENSION(1:3,1:Nstep, 1:Nprobe) :: probe = 0.0_dp

    CALL set_up_sim(material, NCeldas, Nx, Ny, Nz)
    CALL allocate_memory() !Mechanical variables
    CALL load_material(material, rho_inv, c_E, beta_s, e_piezo)
    CALL load_PML() !w1, w2
    CALL open_vtk_file(outfile)
    CALL write_point_data(outfile, data_name)


    WRITE(*, *) w1(UROLL3(2, 2, 3), z)

    DO step = 1, Nstep
        CALL v_half_step()
        CALL free_boundary_v()
        CALL dot_source()
        CALL T_half_step()
        CALL free_boundary_T()
        IF (MOD(STEP, 100) .EQ. 0) THEN
            1000 format('onda', i3.3, '_'i3.3, '.vtk')
            WRITE(outfile, 1000) 1, step/100
            CALL open_vtk_file(outfile)
            CALL write_point_data(outfile, data_name)
        END IF
        cont = 1
        DO ix = Nx/2-16, Nx/2+16, 4
            DO iy = Ny/2-16, Ny/2+16, 4
                DO iz = 1, 2*8+1, 2
                    probe(:, step, cont) = (/ Vx(UROLL3(ix, iy, iz)), Vy(UROLL3(ix, iy, iz)), Vz(UROLL3(ix, iy, iz)) /)
                    cont = cont + 1
                END DO
            END DO
        END DO
        write(*, *) step
    END DO

!~     OPEN(UNIT = 12, FILE = 'probe', ACTION = "write", STATUS = "replace", form = 'unformatted')
!~     write(12) probe
!~     CLOSE(12)

    DEALLOCATE(c_E)
    DEALLOCATE(beta_s)
    DEALLOCATE(e_piezo)

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

    DEALLOCATE(w1)
    DEALLOCATE(w2)


END PROGRAM
