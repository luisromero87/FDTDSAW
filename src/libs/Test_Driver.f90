Program Test_LIB_VTK_IO
USE IR_Precision  
USE LIB_VTK_IO
implicit none

! variabili per l'esempio di "RectilinearGrid" XML
integer(I4P), parameter:: nx1=0_I4P
integer(I4P), parameter:: nx2=1_I4P
integer(I4P), parameter:: ny1=0_I4P
integer(I4P), parameter:: ny2=1_I4P
integer(I4P), parameter:: nz1=0_I4P
integer(I4P), parameter:: nz2=1_I4P
real(R8P),    dimension(nx1:nx2):: x_xml_rect
real(R8P),    dimension(ny1:ny2):: y_xml_rect
real(R8P),    dimension(nz1:nz2):: z_xml_rect
integer(I4P), dimension(1:(nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1)):: var_xml_rect_point

! variabili per l'esempio di "STRUCTURED_POINTS"
integer(I4P),   parameter::             Nx   = 3_I4P
integer(I4P),   parameter::             Ny   = 4_I4P
integer(I4P),   parameter::             Nz   = 6_I4P
real(R8P),      parameter::             Dx   = 1._R8P
real(R8P),      parameter::             Dy   = 1._R8P
real(R8P),      parameter::             Dz   = 1._R8P
real(R8P),      parameter::             X0   = 0._R8P
real(R8P),      parameter::             Y0   = 0._R8P
real(R8P),      parameter::             Z0   = 0._R8P
integer(I4P),   dimension(1:Nx*Ny*Nz):: var_str_point

! variabili per l'esempio di "UNSTRUCTURED_GRID"
integer(I4P), parameter::       Nn   = 27_I4P
integer(I4P), parameter::       Ne   = 11_I4P
real(R4P),    dimension(1:Nn):: x_uns
real(R4P),    dimension(1:Nn):: y_uns
real(R4P),    dimension(1:Nn):: z_uns
integer(I4P), dimension(1:Ne):: tipo
integer(I4P), dimension(1:60):: connect
integer(I4P), dimension(1:49):: connect_xml
integer(I4P), dimension(1:Ne):: offset_xml
real(R8P),    dimension(1:Nn):: var_uns_grid
integer(I4P), dimension(1:Nn):: var_uns_grid_X
integer(I4P), dimension(1:Nn):: var_uns_grid_Y
integer(I4P), dimension(1:Nn):: var_uns_grid_Z

integer(I4P):: E_IO

! esempio di output in "RectilinearGrid" XML (binario)
! creazione del file
E_IO = VTK_INI_XML(output_format = 'BINARY',              &
                   filename      = 'XML_RECT_BINARY.vtr', &
                   mesh_topology = 'RectilinearGrid',     &
                   nx1=nx1,nx2=nx2,ny1=ny1,ny2=ny2,nz1=nz1,nz2=nz2)
! salvataggio della geometria del "Piece" corrente
x_xml_rect=(/0._R8P,1._R8P/)
y_xml_rect=(/0._R8P,1._R8P/)
z_xml_rect=(/0._R8P,1._R8P/)
E_IO = VTK_GEO_XML(nx1=nx1,nx2=nx2,ny1=ny1,ny2=ny2,nz1=nz1,nz2=nz2, &
	                 X=x_xml_rect,Y=y_xml_rect,Z=z_xml_rect)
! inizializzazione del salvataggio delle variabili definite ai nodi del "Piece" corrente
E_IO = VTK_DAT_XML(var_location     = 'node', &
                   var_block_action = 'OPEN')
! salvataggio delle variabili definite ai nodi del "Piece" corrente
var_xml_rect_point(1:2)=(/0_I4P,1_I4P/)
var_xml_rect_point(3:4)=(/7_I4P,5_I4P/)
var_xml_rect_point(5:6)=(/2_I4P,6_I4P/)
var_xml_rect_point(7:8)=(/9_I4P,0_I4P/)
E_IO = VTK_VAR_XML(NC_NN   = (nx2-nx1+1)*(ny2-ny1+1)*(nz2-nz1+1), &
                   varname = 'volume_scalars',                    &
                   var     = var_xml_rect_point)
! chiusura del blocco delle variabili definite ai nodi del "Piece" corrente
E_IO = VTK_DAT_XML(var_location     = 'node', &
                   var_block_action = 'Close')
! chiusura del "Piece" corrente
E_IO = VTK_GEO_XML()
! chiusura del file
E_IO = VTK_END_XML()

! esempio di output in STRUCTURED_POINTS (ascii)
! creazione del file
E_IO = VTK_INI(output_format = 'ASCII',                     &
               filename      = 'STR_POINTS_ASCII.vtk',      &
               title         = 'Structured Points Example', &
               mesh_topology = 'STRUCTURED_POINTS')
! salvataggio della geometria
E_IO = VTK_GEO(Nx=Nx,Ny=Ny,Nz=Nz, &
               X0=X0,Y0=Y0,Z0=Z0, &
               Dx=Dx,Dy=Dy,Dz=Dz)
! inizializzazione del salvataggio delle variabili definite ai nodi
E_IO = VTK_DAT(NC_NN        = Nx*Ny*Nz, &
               var_location = 'node')
! salvataggio delle variabili definite ai nodi
var_str_point(1:12)  = (/0,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0/)
var_str_point(13:24) = (/0,5 ,10,15,20,25,25,20,15,10,5 ,0/)
var_str_point(25:36) = (/0,10,20,30,40,50,50,40,30,20,10,0/)
var_str_point(37:48) = (/0,10,20,30,40,50,50,40,30,20,10,0/)
var_str_point(49:60) = (/0,5 ,10,15,20,25,25,20,15,10,5 ,0/)
var_str_point(61:72) = (/0,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0/)
E_IO = VTK_VAR(NC_NN   = Nx*Ny*Nz,         &
               varname = 'volume_scalars', &
               var     = var_str_point)
! chiusura del file
E_IO = VTK_END()

! esempio di output in UNSTRUCTURED_GRID binario
! creazione del file
E_IO = VTK_INI(output_format = 'BINARY',                    &
               filename      = 'UNST_GRID_BIN.vtk',         &
               title         = 'Unstructured Grid Example', &
               mesh_topology = 'UNSTRUCTURED_GRID')
! salvataggio della geometria
x_uns=(/0,1,2,0,1,2, &
        0,1,2,0,1,2, &
        0,1,2,0,1,2, &
        0,1,2,0,1,2, &
        0,1,2/)
y_uns=(/0,0,0,1,1,1, &
        0,0,0,1,1,1, &
        1,1,1,1,1,1, &
        1,1,1,1,1,1, &
        1,1,1/)
z_uns=(/0,0,0,0,0,0, &
        1,1,1,1,1,1, &
        2,2,2,3,3,3, &
        4,4,4,5,5,5, &
        6,6,6/)
E_IO = VTK_GEO(NN = Nn, &
               X=x_uns,Y=y_uns,Z=z_uns)
! salvataggio della connettività
 connect = (/8,0 ,1 ,4 ,3 ,6 ,7 ,10,9 , &
             8,1 ,2 ,5 ,4 ,7 ,8 ,11,10, &
             4,6 ,10,9 ,12,             &
             4,5 ,11,10,14,             &
             6,15,16,17,14,13,12,       &
             6,18,15,19,16,20,17,       &
             4,22,23,20,19,             &
             3,21,22,18,                &
             3,22,19,18,                &
             2,26,25,                   &
             1,24/)
tipo = (/12, &
         12, &
         10, &
         10, &
         7,  &
         6,  &
         9,  &
         5,  &
         5,  &
         3,  &
         1/)
E_IO = VTK_CON(NC        = Ne,      &
               connect   = connect, &
               cell_type = tipo)
! inizializzazione del salvataggio delle variabili definite ai nodi
E_IO = VTK_DAT(NC_NN        = Nn, &
               var_location = 'node')
! salvataggio delle variabili definite ai nodi
var_uns_grid =(/0.0 ,1.0 ,2.0 ,3.0 ,4.0 ,5.0 , &
                6.0 ,7.0 ,8.0 ,9.0 ,10.0,11.0, &
                12.0,13.0,14.0,15.0,16.0,17.0, &
                18.0,19.0,20.0,21.0,22.0,23.0, &
                24.0,25.0,26.0/)
! variabile scalare
E_IO = VTK_VAR(NC_NN   = Nn, &
               varname = 'scalars', &
               var     = var_uns_grid)
! variabile vettoriale
var_uns_grid_X=(/1,1,0,1,1,0, &
                 1,1,0,1,1,0, &
                 0,0,0,0,0,0, &
                 0,0,0,0,0,0, &
                 0,0,0/)
var_uns_grid_Y=(/0,1,2,0,1,2, &
                 0,1,2,0,1,2, &
                 0,0,0,0,0,0, &
                 0,0,0,0,0,0, &
                 0,0,0/)
var_uns_grid_Z=(/0,0,0,0,0,0, &
                 0,0,0,0,0,0, &
                 1,1,1,1,1,1, &
                 1,1,1,1,1,1, &
                 1,1,1/)
E_IO = VTK_VAR(NC_NN   = Nn,        &
               varname = 'vectors', &
               varX=var_uns_grid_X,varY=var_uns_grid_Y,varZ=var_uns_grid_Z)
! chiusura del file
E_IO = VTK_END()

! esempio di output in "UnstructuredGrid" XML (binario)
! creazione del file
E_IO = VTK_INI_XML(output_format = 'BINARY',              &
                   filename      = 'XML_UNST_BINARY.vtu', &
                   mesh_topology = 'UnstructuredGrid')
! salvataggio della geometria del "Piece" corrente
E_IO = VTK_GEO_XML(NN = Nn, &
	                 NC = Ne, &
	                 X=x_uns,Y=y_uns,Z=z_uns)
! salvataggio della connettività
 connect_xml = (/0 ,1 ,4 ,3 ,6 ,7 ,10,9 , &
                 1 ,2 ,5 ,4 ,7 ,8 ,11,10, &
                 6 ,10,9 ,12,             &
                 5 ,11,10,14,             &
                 15,16,17,14,13,12,       &
                 18,15,19,16,20,17,       &
                 22,23,20,19,             &
                 21,22,18,                &
                 22,19,18,                &
                 26,25,                   &
                 24/)
offset_xml = (/8 , &
               16, &
               20, &
               24, &
               30, &
               36, &
               40, &
               43, &
               46, &
               48, &
               49/)
E_IO = VTK_CON_XML(NC        = Ne,          &
                   connect   = connect_xml, &
                   offset    = offset_xml,  &
                   cell_type = (/12_I1P,    &
                                 12_I1P,    &
                                 10_I1P,    &
                                 10_I1P,    &
                                 7_I1P,     &
                                 6_I1P,     &
                                 9_I1P,     &
                                 5_I1P,     &
                                 5_I1P,     &
                                 3_I1P,     &
                                 1_I1P/))
! inizializzazione del salvataggio delle variabili definite ai nodi del "Piece" corrente
E_IO = VTK_DAT_XML(var_location     = 'node', &
                   var_block_action = 'opeN')
! salvataggio delle variabili definite ai nodi del "Piece" corrente
! variabile scalare
E_IO = VTK_VAR_XML(NC_NN   = Nn,        &
                   varname = 'scalars', &
                   var     = var_uns_grid)
! variabile vettoriale
E_IO = VTK_VAR_XML(NC_NN   = Nn,       &
                   varname = 'vector', &
                   varX=var_uns_grid_X,varY=var_uns_grid_Y,varZ=var_uns_grid_Z)
! chiusura del blocco delle variabili definite ai nodi del "Piece" corrente
E_IO = VTK_DAT_XML(var_location     = 'node', &
                   var_block_action = 'CLOSE')
! chiusura del "Piece" corrente
E_IO = VTK_GEO_XML()
! chiusura del file
E_IO = VTK_END_XML()


! secondo esempio di output in "UnstructuredGrid" XML (binario)
! creazione del file
E_IO = VTK_INI_XML(output_format = 'BINARY',                &
                   filename      = 'XML_UNST_BINARY-2.vtu', &
                   mesh_topology = 'UnstructuredGrid')
! salvataggio della geometria del "Piece" corrente
E_IO = VTK_GEO_XML(NN = 4,                            &
	                 NC = 1,                            &
	                 X=(/0._R8P,1._R8P,1._R8P,0._R8P/), &
	                 Y=(/0._R8P,0._R8P,1._R8P,1._R8P/), &
	                 Z=(/0._R8P,0._R8P,0._R8P,0._R8P/))
! salvataggio della connettività
connect_xml(1:4) = (/0, 1, 2, 3/)
E_IO = VTK_CON_XML(NC        = 1,                &
                   connect   = connect_xml(1:4), &
                   offset    = (/4/),            &
                   cell_type = (/9_I1P/))
! inizializzazione del salvataggio delle variabili definite ai nodi del "Piece" corrente
E_IO = VTK_DAT_XML(var_location     = 'node', &
                   var_block_action = 'OPEN')
! salvataggio delle variabili definite ai nodi del "Piece" corrente
! variabile scalare
E_IO = VTK_VAR_XML(NC_NN   = 4,         &
                   varname = 'scalars', &
                   var     = (/1,2,3,4/))
! chiusura del blocco delle variabili definite ai nodi del "Piece" corrente
E_IO = VTK_DAT_XML(var_location     = 'node', &
                   var_block_action = 'CLOSE')
! chiusura del "Piece" corrente
E_IO = VTK_GEO_XML()
! chiusura del file
E_IO = VTK_END_XML()

endprogram Test_LIB_VTK_IO
