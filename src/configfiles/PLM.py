from numpy import *
import struct


def find_parameter(parameter):
    f = open('sim.config', 'r')
    buff=''
    while buff.find(parameter):
        buff=f.readline()
        buff.find(parameter)

    for subbuff in buff.split('\t'):
	    try:
	        value= float(subbuff)
	    except:
		    pass
		    
    f.close()
    return value

def fomega(i,smax,PLM_w,m):
	return smax*((i+1.0)/PLM_w)**m
	
Nx=int(find_parameter('Nx'))
Ny=int(find_parameter('Ny'))
Nz=int(find_parameter('Nz'))
dt=find_parameter('dt')
dx=find_parameter('dx')
dy=find_parameter('dy')
dz=find_parameter('dz')

Omega_x=zeros((Nx,Ny,Nz),dtype='float64')
Omega_y=zeros((Nx,Ny,Nz),dtype='float64')
Omega_z=zeros((Nx,Ny,Nz),dtype='float64')
m=1.55
PLM_w=10
smax=0.0041*PLM_w/dt

#PLM  <- x ->
for ix in range(PLM_w):
    for iy in range(Ny):
        for iz in range(Nz):
			omega=fomega(ix,smax,PLM_w,m)
			#<- x
			Omega_x[PLM_w-1-ix,iy,iz]=omega
            #Omega_y[PLM_w-1-ix,iy,iz]=smax*(ix**2.0)
            #Omega_z[PLM_w-1-ix,iy,iz]=smax*(ix**2.0)
			#   x ->
			Omega_x[Nx-PLM_w+ix,iy,iz]=omega
            #Omega_y[Nx-PLM_w+ix,iy,iz]=smax*(ix**2.0)
            #Omega_z[Nx-PLM_w+ix,iy,iz]=smax*(ix**2.0)


#PLM  <- y ->
for ix in range(Nx):
    for iy in range(PLM_w):
        for iz in range(Nz):
            omega=fomega(iy,smax,PLM_w,m)
			#<- y
            #Omega_x[ix,PLM_w-1-iy,iz]=sqrt(Omega_x[ix,PLM_w-1-iy,iz]**2+(smax*((iy**2.0)))**2)
            Omega_y[ix,PLM_w-1-iy,iz]=sqrt(Omega_y[ix,PLM_w-1-iy,iz]**2+(omega)**2)
            #Omega_z[ix,PLM_w-1-iy,iz]=sqrt(Omega_z[ix,PLM_w-1-iy,iz]**2+(smax*((iy**2.0)))**2)
            #   y ->
            #Omega_x[ix,Ny-PLM_w+iy,iz]=sqrt(Omega_x[ix,Ny-PLM_w+iy,iz]**2+(smax*((iy**2.0)))**2)
            Omega_y[ix,Ny-PLM_w+iy,iz]=sqrt(Omega_y[ix,Ny-PLM_w+iy,iz]**2+(omega)**2)
            #Omega_z[ix,Ny-PLM_w+iy,iz]=sqrt(Omega_z[ix,Ny-PLM_w+iy,iz]**2+(smax*((iy**2.0)))**2)
	  
#PLM     z ->
for ix in range(Nx):
    for iy in range(Ny):
        for iz in range(PLM_w):
            omega=fomega(iz,smax,PLM_w,m)
			#   y ->
			#Omega_x[ix,iy,Nz-PLM_w+iz]=sqrt(Omega_x[ix,iy,Nz-PLM_w+iz]**2+(smax*((iz**2.0)))**2)
			#Omega_y[ix,iy,Nz-PLM_w+iz]=sqrt(Omega_y[ix,iy,Nz-PLM_w+iz]**2+(smax*((iz**2.0)))**2)
            Omega_z[ix,iy,Nz-PLM_w+iz]=sqrt(Omega_z[ix,iy,Nz-PLM_w+iz]**2+(omega)**2)

w1_x=zeros((Nx,Ny,Nz),dtype='float64')
w1_y=zeros((Nx,Ny,Nz),dtype='float64')
w1_z=zeros((Nx,Ny,Nz),dtype='float64')
w2_x=zeros((Nx,Ny,Nz),dtype='float64')
w2_y=zeros((Nx,Ny,Nz),dtype='float64')
w2_z=zeros((Nx,Ny,Nz),dtype='float64')

for ix in range(Nx):
    for iy in range(Ny):
        for iz in range(Nz):
            w1_x[ix,iy,iz] = (2.0-dt*Omega_x[ix,iy,iz])/(2.0+dt*Omega_x[ix,iy,iz])
            w1_y[ix,iy,iz] = (2.0-dt*Omega_y[ix,iy,iz])/(2.0+dt*Omega_y[ix,iy,iz])
            w1_z[ix,iy,iz] = (2.0-dt*Omega_z[ix,iy,iz])/(2.0+dt*Omega_z[ix,iy,iz])
            w2_x[ix,iy,iz] = 2.0*dt /( dx*(2.0+dt*Omega_x[ix,iy,iz]) )
            w2_y[ix,iy,iz] = 2.0*dt /( dy*(2.0+dt*Omega_y[ix,iy,iz]) )
            w2_z[ix,iy,iz] = 2.0*dt /( dz*(2.0+dt*Omega_z[ix,iy,iz]) )


f=open('PLM','wb') 

rlc=8*Ny*Nz

for row in w1_x[:]:
    f.write(struct.pack('i',rlc))
    (row.transpose()).tofile(f)
    f.write(struct.pack('i',rlc))


for row in w1_y[:]:
    f.write(struct.pack('i',rlc))
    (row.transpose()).tofile(f)
    f.write(struct.pack('i',rlc))


for row in w1_z[:]:
    f.write(struct.pack('i',rlc))
    (row.transpose()).tofile(f)
    f.write(struct.pack('i',rlc))



for row in w2_x[:]:
    f.write(struct.pack('i',rlc))
    (row.transpose()).tofile(f)
    f.write(struct.pack('i',rlc))


for row in w2_y[:]:
    f.write(struct.pack('i',rlc))
    (row.transpose()).tofile(f)
    f.write(struct.pack('i',rlc))


for row in w2_z[:]:
    f.write(struct.pack('i',rlc))
    (row.transpose()).tofile(f)
    f.write(struct.pack('i',rlc))

f.close()





