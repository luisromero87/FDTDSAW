# Environment
MKDIR=mkdir
CP=cp
RM=rm
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=cc
CCC=CC
CXX=CC
FC=f95
AS=as

#Compiler
FC=gfortran
#Flags
FFLAGS= -O3

# Macros
CND_PLATFORM=Linux-x86
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build
CND_SRCDIR=src
CND_LIBSDIR=libs

#Source Directory
SRCDIR=${CND_SRCDIR}
SRCLIBSDIR=${CND_SRCDIR}/${CND_LIBSDIR}

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/fdtdsaw.o \
	${OBJECTDIR}/libs/type_kinds.o \
	${OBJECTDIR}/libs/constants_module.o \
	${OBJECTDIR}/libs/global_vars.o \
	${OBJECTDIR}/libs/set_up_sim.o \
	${OBJECTDIR}/libs/uroll3.o \
	${OBJECTDIR}/libs/allocate_memory.o \
	${OBJECTDIR}/libs/load_material.o \
	${OBJECTDIR}/libs/load_pml.o \
	${OBJECTDIR}/libs/write_to_vtk.o \
	${OBJECTDIR}/libs/time_step.o \

#BUILD MAIN
.PHONY: FDTDSAW
FDTDSAW: ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fdtdsaw
${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fdtdsaw: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	$(FC) -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fdtdsaw ${OBJECTFILES} ${LDLIBSOPTIONS} 

#BUILD MAIN
.PHONY: main
main: ${OBJECTDIR}/fdtdsaw.o
${OBJECTDIR}/fdtdsaw.o: ${SRCDIR}/main.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(FC) -c -o ${OBJECTDIR}/fdtdsaw.o ${SRCDIR}/main.f90
	
#BUILD LIBS
.PHONY: type_kinds
type_kinds: ${OBJECTDIR}/libs/type_kinds.o
${OBJECTDIR}/libs/type_kinds.o: ${SRCLIBSDIR}/type_kinds.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) -c -o ${OBJECTDIR}/libs/type_kinds.o ${SRCLIBSDIR}/type_kinds.f90

.PHONY: constants_module
constants_module: ${OBJECTDIR}/libs/constants_module.o \
                  type_kinds
${OBJECTDIR}/libs/constants_module.o: ${SRCLIBSDIR}/constants_module.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) -c -o ${OBJECTDIR}/libs/constants_module.o ${SRCLIBSDIR}/constants_module.f90

.PHONY: global_vars
global_vars: ${OBJECTDIR}/libs/global_vars.o \
                  type_kinds \
                  constants_module
${OBJECTDIR}/libs/global_vars.o: ${SRCLIBSDIR}/global_vars.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) -c -o ${OBJECTDIR}/libs/global_vars.o ${SRCLIBSDIR}/global_vars.f90

.PHONY: set_up_sim
set_up_sim: ${OBJECTDIR}/libs/set_up_sim.o \
                  type_kinds \
                  constants_module
${OBJECTDIR}/libs/set_up_sim.o: ${SRCLIBSDIR}/set_up_sim.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) -c -o ${OBJECTDIR}/libs/set_up_sim.o ${SRCLIBSDIR}/set_up_sim.f90

.PHONY: uroll3
uroll3: ${OBJECTDIR}/libs/uroll3.o \
                  type_kinds \
                  constants_module \
                  global_vars
${OBJECTDIR}/libs/uroll3.o: ${SRCLIBSDIR}/uroll3.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) -c -o ${OBJECTDIR}/libs/uroll3.o ${SRCLIBSDIR}/uroll3.f90

.PHONY: allocate_memory
allocate_memory: ${OBJECTDIR}/libs/allocate_memory.o \
                  type_kinds \
                  constants_module \
                  global_vars
${OBJECTDIR}/libs/allocate_memory.o: ${SRCLIBSDIR}/allocate_memory.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) -c -o ${OBJECTDIR}/libs/allocate_memory.o ${SRCLIBSDIR}/allocate_memory.f90

.PHONY: load_material
load_material: ${OBJECTDIR}/libs/load_material.o \
                  type_kinds \
                  constants_module \
                  global_vars
${OBJECTDIR}/libs/load_material.o: ${SRCLIBSDIR}/load_material.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) -c -o ${OBJECTDIR}/libs/load_material.o ${SRCLIBSDIR}/load_material.f90

.PHONY: load_pml
load_pml: ${OBJECTDIR}/libs/load_pml.o \
                  type_kinds \
                  constants_module \
                  global_vars
${OBJECTDIR}/libs/load_pml.o: ${SRCLIBSDIR}/load_pml.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) -c -o ${OBJECTDIR}/libs/load_pml.o ${SRCLIBSDIR}/load_pml.f90

.PHONY: write_to_vtk
write_to_vtk: ${OBJECTDIR}/libs/write_to_vtk.o \
                  type_kinds \
                  constants_module \
                  global_vars
${OBJECTDIR}/libs/write_to_vtk.o: ${SRCLIBSDIR}/write_to_vtk.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) -c -o ${OBJECTDIR}/libs/write_to_vtk.o ${SRCLIBSDIR}/write_to_vtk.f90

.PHONY: time_step
time_step: ${OBJECTDIR}/libs/time_step.o \
                  type_kinds \
                  constants_module \
                  global_vars
${OBJECTDIR}/libs/time_step.o: ${SRCLIBSDIR}/time_step.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) -c -o ${OBJECTDIR}/libs/time_step.o ${SRCLIBSDIR}/time_step.f90

#CLEAN
.PHONY: clean
clean:
	${RM} -r ${OBJECTDIR}


#~ FDTDSAW: $(OBJECTS)
	#~ $(FC) -o FDTDSAW $(OBJECTS) $(FFLAGS)
#~ 
#~ libs/%.o: %.f90
	#~ $(FC) $(FFLAGS) -I/libs -L/libs -c $<


# include project implementation makefile
#include nbproject/Makefile-impl.mk

# include project make variables
#include nbproject/Makefile-variables.mk
