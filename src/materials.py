from numpy import *
import struct

#---------------------
#      Bi12GeO20      |
#---------------------
f=open('Materials/bi12geo20','wb') 

# Density
rho=9200.0
rho_inv=array(1.0/rho,dtype='float64')
rlc=8*1
f.write(struct.pack('i',8))
rho_inv.tofile(f)
f.write(struct.pack('i',8))

# Stiffness at constant E
c_E=zeros((6,6),dtype='float64')

c_E[0,0]=c_E[1,1]=c_E[2,2]=12.8e10
c_E[0,1]=c_E[0,2]=c_E[1,2]=c_E[1,0]=c_E[2,0]=c_E[2,1]=3.05e10
c_E[3,3]=c_E[4,4]=c_E[5,5]=2.55e10

rlc=8*6*6

f.write(struct.pack('i',rlc))
c_E.tofile(f)
f.write(struct.pack('i',rlc))

# Compliance at constant E
s_E=linalg.inv(c_E)

rlc=8*6*6

f.write(struct.pack('i',rlc))
s_E.tofile(f)
f.write(struct.pack('i',rlc))

#Electric permitivity 
eps_s=zeros((3,3),dtype='float64')
beta_s=zeros((3,3),dtype='float64')

eps_s[0,0]=34.2e-11
eps_s[1,1]=34.2e-11
eps_s[2,2]=34.2e-11

beta_s=linalg.inv(eps_s)

rlc=8*3*3

f.write(struct.pack('i',rlc))
beta_s.tofile(f)
f.write(struct.pack('i',rlc))


#Dielectric constant
e_piezo=zeros((3,6),dtype='float64')

e_piezo[0,3]=e_piezo[1,4]=e_piezo[2,5]=0.99

rlc=8*3*6

f.write(struct.pack('i',rlc))
(e_piezo.transpose()).tofile(f)
f.write(struct.pack('i',rlc))

f.close()

s=zeros((6,6))
d=zeros((3,6))


#---------------------
#      Bi4Ge3O12      |
#---------------------
f=open('Materials/bi4ge3o12','wb') 

# Density
rho=7095.0
rho_inv=array(1.0/rho,dtype='float64')
rlc=8*1
f.write(struct.pack('i',8))
rho_inv.tofile(f)
f.write(struct.pack('i',8))

# Stiffness at constant E
c_E=zeros((6,6),dtype='float64')

c_E[0,0]=c_E[1,1]=c_E[2,2]=11.58e10
c_E[0,1]=c_E[0,2]=c_E[1,2]=c_E[1,0]=c_E[2,0]=c_E[2,1]=2.07e10
c_E[3,3]=c_E[4,4]=c_E[5,5]=4.36e10

rlc=8*6*6

f.write(struct.pack('i',rlc))
c_E.tofile(f)
f.write(struct.pack('i',rlc))

# Compliance at constant E
s_E=linalg.inv(c_E)

rlc=8*6*6

f.write(struct.pack('i',rlc))
s_E.tofile(f)
f.write(struct.pack('i',rlc))

#Electric permitivity 
eps_s=zeros((3,3),dtype='float64')
beta_s=zeros((3,3),dtype='float64')

eps_s[0,0]=eps_s[1,1]=eps_s[2,2]=16*8.8541878176e-12

beta_s=linalg.inv(eps_s)

rlc=8*3*3

f.write(struct.pack('i',rlc))
beta_s.tofile(f)
f.write(struct.pack('i',rlc))


#Dielectric constant
e_piezo=zeros((3,6),dtype='float64')

e_piezo[0,3]=e_piezo[1,4]=e_piezo[2,5]=0.0376

rlc=8*3*6

f.write(struct.pack('i',rlc))
(e_piezo.transpose()).tofile(f)
f.write(struct.pack('i',rlc))

f.close()

s=zeros((6,6))
d=zeros((3,6))


#--------------------------
#  Barium Sodium Niobate  |
#--------------------------
f=open('Materials/ba2nanb5o15','wb') 

# Density
rho=5300.0
rho_inv=array(1.0/rho,dtype='float64')
rlc=8*1
f.write(struct.pack('i',8))
rho_inv.tofile(f)
f.write(struct.pack('i',8))

# Stiffness at constant E
c_E=zeros((6,6),dtype='float64')

c_E[0,0] = 23.9e10
c_E[1,1] = 24.7e10
c_E[2,2] = 13.5e10
c_E[0,1] = c_E[1,0] = 10.4e10
c_E[0,2] = c_E[2,0] = 5.0e10
c_E[1,2] = c_E[2,1] = 5.2e10
c_E[3,3] = 6.5e10
c_E[4,4] = 6.6e10
c_E[5,5] = 7.6e10

rlc=8*6*6

f.write(struct.pack('i',rlc))
c_E.tofile(f)
f.write(struct.pack('i',rlc))

# Stiffness at constant E
s_E=zeros((6,6),dtype='float64')

s_E[0,0] = 5.3e-12
s_E[1,1] = 5.14e-12
s_E[2,2] = 8.33e-12
s_E[0,1] = s_E[1,0] = -1.98e-12
s_E[0,2] = s_E[2,0] = -1.2e-12
s_E[1,2] = s_E[2,1] = -1.25e-12
s_E[3,3] = 15.4e-12
s_E[4,4] = 15.2e-12
s_E[5,5] = 13.2e-12

rlc=8*6*6

f.write(struct.pack('i',rlc))
s_E.tofile(f)
f.write(struct.pack('i',rlc))

#Electric permitivity 
eps_s=zeros((3,3),dtype='float64')
beta_s=zeros((3,3),dtype='float64')

eps_s[0,0]=222
eps_s[1,1]=227
eps_s[2,2]=32
eps_s=eps_s*8.8541878176e-12

beta_s=linalg.inv(eps_s)

rlc=8*3*3

f.write(struct.pack('i',rlc))
beta_s.tofile(f)
f.write(struct.pack('i',rlc))


#Piezoelectric stress constants
e_piezo=zeros((3,6),dtype='float64')

e_piezo[0,4]=2.8
e_piezo[1,3]=3.4
e_piezo[2,0]=-0.4
e_piezo[2,1]=-0.3
e_piezo[2,2]=4.3

rlc=8*3*6

f.write(struct.pack('i',rlc))
(e_piezo.transpose()).tofile(f)
f.write(struct.pack('i',rlc))

f.close()


#--------------------------
#       Rutile            |
#--------------------------
f=open('Materials/rutile','wb') 

# Density
rho=4260.0
rho_inv=array(1.0/rho,dtype='float64')
rlc=8*1
f.write(struct.pack('i',8))
rho_inv.tofile(f)
f.write(struct.pack('i',8))

# Stiffness at constant E
c_E=zeros((6,6),dtype='float64')

c_E[0,0] = 26.6e10
c_E[1,1] = 26.6e10
c_E[2,2] = 46.99e10
c_E[0,1] = c_E[1,0] = 17.33e10
c_E[0,2] = c_E[2,0] = 13.62e10
c_E[1,2] = c_E[2,1] = 13.62e10
c_E[3,3] = 12.39e10
c_E[4,4] = 12.39e10
c_E[5,5] = 18.86e10

rlc=8*6*6

f.write(struct.pack('i',rlc))
c_E.tofile(f)
f.write(struct.pack('i',rlc))

# Compliance at constant E
s_E=linalg.inv(c_E)

rlc=8*6*6

f.write(struct.pack('i',rlc))
s_E.tofile(f)
f.write(struct.pack('i',rlc))

#Electric permitivity 
eps_s=zeros((3,3),dtype='float64')
beta_s=zeros((3,3),dtype='float64')

eps_s[0,0]=34.2e-11/38*222*0
eps_s[1,1]=34.2e-11/38*227*0
eps_s[2,2]=34.2e-11/38*32*0

#beta_s=linalg.inv(eps_s)
beta_s=eps_s

rlc=8*3*3

f.write(struct.pack('i',rlc))
beta_s.tofile(f)
f.write(struct.pack('i',rlc))


#Dielectric constant
e_piezo=zeros((3,6),dtype='float64')

e_piezo[0,3]=e_piezo[1,4]=e_piezo[2,5]=0.99*0

rlc=8*3*6

f.write(struct.pack('i',rlc))
(e_piezo.transpose()).tofile(f)
f.write(struct.pack('i',rlc))

f.close()

s=zeros((6,6))
d=zeros((3,6))

#--------------------------
#    Cadmium Sulfide      |
#--------------------------
f=open('Materials/cds','wb') 

# Density
rho=4820.0
rho_inv=array(1.0/rho,dtype='float64')
rlc=8*1
f.write(struct.pack('i',8))
rho_inv.tofile(f)
f.write(struct.pack('i',8))

# Stiffness at constant E
c_E=zeros((6,6),dtype='float64')

c_E[0,0] = 9.07e10
c_E[1,1] = 9.07e10
c_E[2,2] = 9.38e10
c_E[0,1] = c_E[1,0] = 5.81e10
c_E[0,2] = c_E[2,0] = 5.10e10
c_E[1,2] = c_E[2,1] = 5.10e10
c_E[3,3] = 1.504e10
c_E[4,4] = 1.504e10
c_E[5,5] = 0.5*(c_E[0,0]-c_E[0,1])

rlc=8*6*6

f.write(struct.pack('i',rlc))
c_E.tofile(f)
f.write(struct.pack('i',rlc))

# Compliance at constant E
s_E=linalg.inv(c_E)

rlc=8*6*6

f.write(struct.pack('i',rlc))
s_E.tofile(f)
f.write(struct.pack('i',rlc))

#Electric permitivity 
eps_s=zeros((3,3),dtype='float64')
beta_s=zeros((3,3),dtype='float64')

eps_s[0,0]=34.2e-11/38*222*0
eps_s[1,1]=34.2e-11/38*227*0
eps_s[2,2]=34.2e-11/38*32*0

#beta_s=linalg.inv(eps_s)
beta_s=eps_s

rlc=8*3*3

f.write(struct.pack('i',rlc))
beta_s.tofile(f)
f.write(struct.pack('i',rlc))


#Dielectric constant
e_piezo=zeros((3,6),dtype='float64')

e_piezo[0,3]=e_piezo[1,4]=e_piezo[2,5]=0.99*0

rlc=8*3*6

f.write(struct.pack('i',rlc))
(e_piezo.transpose()).tofile(f)
f.write(struct.pack('i',rlc))

f.close()

s=zeros((6,6))
d=zeros((3,6))
