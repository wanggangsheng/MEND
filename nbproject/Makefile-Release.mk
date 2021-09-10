#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-MacOSX
CND_DLIB_EXT=dylib
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/MEND_IN.o \
	${OBJECTDIR}/src/MEND_main.o \
	${OBJECTDIR}/src/MOD_MCMC.o \
	${OBJECTDIR}/src/MOD_MEND.o \
	${OBJECTDIR}/src/MOD_MEND_TYPE.o \
	${OBJECTDIR}/src/MOD_OPT.o \
	${OBJECTDIR}/src/MOD_OPT_TYPE.o \
	${OBJECTDIR}/src/MOD_STRING.o \
	${OBJECTDIR}/src/MOD_USRFS.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mend

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mend: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.f} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mend ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/src/MEND_IN.o: src/MEND_IN.F90
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/MEND_IN.o src/MEND_IN.F90

${OBJECTDIR}/src/MEND_main.o: src/MEND_main.F90
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/MEND_main.o src/MEND_main.F90

${OBJECTDIR}/src/MOD_MCMC.o: src/MOD_MCMC.F90
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/MOD_MCMC.o src/MOD_MCMC.F90

${OBJECTDIR}/src/MOD_MEND.o: src/MOD_MEND.F90
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/MOD_MEND.o src/MOD_MEND.F90

${OBJECTDIR}/src/MOD_MEND_TYPE.o: src/MOD_MEND_TYPE.F90
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/MOD_MEND_TYPE.o src/MOD_MEND_TYPE.F90

${OBJECTDIR}/src/MOD_OPT.o: src/MOD_OPT.F90
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/MOD_OPT.o src/MOD_OPT.F90

${OBJECTDIR}/src/MOD_OPT_TYPE.o: src/MOD_OPT_TYPE.F90
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/MOD_OPT_TYPE.o src/MOD_OPT_TYPE.F90

${OBJECTDIR}/src/MOD_STRING.o: src/MOD_STRING.F90
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/MOD_STRING.o src/MOD_STRING.F90

${OBJECTDIR}/src/MOD_USRFS.o: src/MOD_USRFS.F90
	${MKDIR} -p ${OBJECTDIR}/src
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/src/MOD_USRFS.o src/MOD_USRFS.F90

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} *.mod

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
