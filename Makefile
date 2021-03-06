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
FC=gfortran
FFLAGS= -O3 -fbackslash -cpp
AS=as
MPI= no
PREPROC=
DEBUG=no

# Macros
CND_PLATFORM=Linux-x86
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build
CND_SRCDIR=src
CND_LIBSDIR=libs
CND_MATERIALSDIR=Materials
CND_CONFIGFILESDIR=configfiles

#Source Directory
SRCDIR=${CND_SRCDIR}
SRCLIBSDIR=${CND_SRCDIR}/${CND_LIBSDIR}
SRCMATERIALSDIR=${CND_SRCDIR}/${CND_MATERIALSDIR}
SRCCONFIGFILESDIR=${CND_SRCDIR}/${CND_CONFIGFILESDIR}

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/libs/IR_Precision.o \
	${OBJECTDIR}/libs/Lib_Base64.o \
	${OBJECTDIR}/libs/Lib_VTK_IO.o \
	${OBJECTDIR}/libs/type_kinds.o \
	${OBJECTDIR}/libs/constants_module.o \
	${OBJECTDIR}/libs/global_vars.o \
	${OBJECTDIR}/libs/time_step.o \
	${OBJECTDIR}/libs/Lib_FDTD_SAW.o \
	${OBJECTDIR}/fdtdsaw.o

ifeq "$(DEBUG)" "yes"
  FFLAGS := $(FFLAGS) -Wall -fcheck=all -C -g -fbacktrace
endif
ifeq "$(MPI)" "yes"
  PREPROC := $(PREPROC) -DMPI2
  FC = mpif90
endif

FFLAGS := $(FFLAGS) $(PREPROC)

#BUILD MAIN
.PHONY: FDTDSAW
FDTDSAW: ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fdtdsaw
${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fdtdsaw: ${OBJECTFILES} 
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	$(FC) $(FFLAGS) -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fdtdsaw ${OBJECTFILES} ${LDLIBSOPTIONS} 
	$(CP) -nr ${SRCMATERIALSDIR} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	$(CP) -ur ${SRCCONFIGFILESDIR} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	$(CP) -nr ${SRCDIR}/timer.sh ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/

#BUILD MAIN
.PHONY: main
main: ${OBJECTDIR}/fdtdsaw.o
${OBJECTDIR}/fdtdsaw.o: ${SRCDIR}/main.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(FC) $(FFLAGS) -c -o ${OBJECTDIR}/fdtdsaw.o ${SRCDIR}/main.f90
	
#BUILD LIBS
.PHONY: type_kinds
type_kinds: ${OBJECTDIR}/libs/type_kinds.o
${OBJECTDIR}/libs/type_kinds.o: ${SRCLIBSDIR}/type_kinds.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) $(FFLAGS) -c -o ${OBJECTDIR}/libs/type_kinds.o ${SRCLIBSDIR}/type_kinds.f90

.PHONY: constants_module
constants_module: type_kinds \
                  ${OBJECTDIR}/libs/constants_module.o
${OBJECTDIR}/libs/constants_module.o: ${SRCLIBSDIR}/constants_module.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) $(FFLAGS) -c -o ${OBJECTDIR}/libs/constants_module.o ${SRCLIBSDIR}/constants_module.f90

.PHONY: global_vars
global_vars:      type_kinds \
                  constants_module \
                  ${OBJECTDIR}/libs/global_vars.o 
${OBJECTDIR}/libs/global_vars.o: ${SRCLIBSDIR}/global_vars.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) $(FFLAGS) -c -o ${OBJECTDIR}/libs/global_vars.o ${SRCLIBSDIR}/global_vars.f90

.PHONY: time_step
time_step:        type_kinds \
                  constants_module \
                  global_vars \
                  ${OBJECTDIR}/libs/time_step.o
${OBJECTDIR}/libs/time_step.o: ${SRCLIBSDIR}/time_step.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) $(FFLAGS) -c -o ${OBJECTDIR}/libs/time_step.o ${SRCLIBSDIR}/time_step.f90
	
.PHONY: Lib_FDTD_SAW
Lib_FDTD_SAW:     Lib_VTK_IO \
                  type_kinds \
                  constants_module \
                  global_vars \
                  ${OBJECTDIR}/libs/Lib_FDTD_SAW.o 
${OBJECTDIR}/libs/Lib_FDTD_SAW.o: ${SRCLIBSDIR}/Lib_FDTD_SAW.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) $(FFLAGS) -c -o ${OBJECTDIR}/libs/Lib_FDTD_SAW.o ${SRCLIBSDIR}/Lib_FDTD_SAW.f90


.PHONY: Lib_VTK_IO
Lib_VTK_IO:       IR_Precision \
                  Lib_Base64 \
                  ${OBJECTDIR}/libs/Lib_VTK_IO.o
${OBJECTDIR}/libs/Lib_VTK_IO.o: ${SRCLIBSDIR}/Lib_VTK_IO.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) $(FFLAGS) -c -o ${OBJECTDIR}/libs/Lib_VTK_IO.o ${SRCLIBSDIR}/Lib_VTK_IO.f90

.PHONY: Lib_Base64
Lib_Base64: IR_Precision \
            ${OBJECTDIR}/libs/Lib_Base64.o
${OBJECTDIR}/libs/Lib_Base64.o: ${SRCLIBSDIR}/Lib_Base64.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) $(FFLAGS) -c -o ${OBJECTDIR}/libs/Lib_Base64.o ${SRCLIBSDIR}/Lib_Base64.f90

.PHONY: IR_Precision
IR_Precision: ${OBJECTDIR}/libs/IR_Precision.o 
${OBJECTDIR}/libs/IR_Precision.o: ${SRCLIBSDIR}/IR_Precision.f90 
	${MKDIR} -p ${OBJECTDIR}/libs
	$(FC) $(FFLAGS) -c -o ${OBJECTDIR}/libs/IR_Precision.o ${SRCLIBSDIR}/IR_Precision.f90
	


#CLEAN
.PHONY: clean
clean:
	${RM} -rf ${OBJECTDIR}
	${RM} -f *.mod
	find ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/ -type d \
	-not -name ${CND_PLATFORM} \
	-not -name 'Materials' -not -name 'configfiles' | xargs rm -rf
