from numpy import *
import struct

#---------------------
#      Bi12GeO20      |
#---------------------
f=open('bi12geo20','wb') 

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

c_E[0,0]=c_E[1,1]=c_E[2,2]=12.8e10
c_E[0,1]=c_E[0,2]=c_E[1,2]=c_E[1,0]=c_E[2,0]=c_E[2,1]=3.05e10
c_E[3,3]=c_E[4,4]=c_E[5,5]=2.55e10

rlc=8*6*6

f.write(struct.pack('i',rlc))
c_E.tofile(f)
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



s=zeros((6,6))
d=zeros((3,6))

