###########################################################################
## Makefile generated for MATLAB file/project 'Calc_pressao'. 
## 
## Makefile     : Calc_pressao_rtw.mk
## Generated on : Mon Aug 08 21:42:03 2022
## MATLAB Coder version: 3.3 (R2017a)
## 
## Build Info:
## 
## Final product: $(RELATIVE_PATH_TO_ANCHOR)/Calc_pressao.lib
## Product type : static-library
## 
###########################################################################

###########################################################################
## MACROS
###########################################################################

# Macro Descriptions:
# PRODUCT_NAME            Name of the system to build
# MAKEFILE                Name of this makefile
# COMPUTER                Computer type. See the MATLAB "computer" command.
# PERL                    PERL Tool
# GEN_LNK_SCRIPT          Perl script to generate the command file
# CMD_FILE                Command file

PRODUCT_NAME              = Calc_pressao
MAKEFILE                  = Calc_pressao_rtw.mk
COMPUTER                  = PCWIN64
MATLAB_ROOT               = C:/PROGRA~1/MATLAB/R2017a
MATLAB_BIN                = C:/PROGRA~1/MATLAB/R2017a/bin
MATLAB_ARCH_BIN           = C:/PROGRA~1/MATLAB/R2017a/bin/win64
MASTER_ANCHOR_DIR         = 
START_DIR                 = D:/Usuarios/IC/C�digos de campo de press�es
ARCH                      = win64
RELATIVE_PATH_TO_ANCHOR   = .
PERL                      = $(MATLAB_ROOT)/sys/perl/win32/bin/perl.exe
GEN_LNK_SCRIPT            = $(MATLAB_ROOT)/rtw/c/tools/mkvc_lnk.pl
CMD_FILE                  = $(PRODUCT_NAME).lnk
C_STANDARD_OPTS           = 
CPP_STANDARD_OPTS         = 

###########################################################################
## TOOLCHAIN SPECIFICATIONS
###########################################################################

# Toolchain Name:          LCC-win64 v2.4.1 | gmake (64-bit Windows)
# Supported Version(s):    2.4.1
# ToolchainInfo Version:   R2017a
# Specification Revision:  1.0
# 

#-----------
# MACROS
#-----------

SHELL              = cmd
LCC_ROOT           = $(MATLAB_ROOT)/sys/lcc64/lcc64
LCC_BUILDLIB       = $(LCC_ROOT)/bin/buildlib
LCC_LIB            = $(LCC_ROOT)/lib64
MW_EXTERNLIB_DIR   = $(MATLAB_ROOT)/extern/lib/win64/microsoft
MW_LIB_DIR         = $(MATLAB_ROOT)/lib/win64
TOOLCHAIN_INCLUDES = -I$(LCC_ROOT)/include64
MEX_OPTS_FILE      = $(MATLAB_ROOT/rtw/c/tools/lcc-win64.xml

TOOLCHAIN_SRCS = 
TOOLCHAIN_INCS = 
TOOLCHAIN_LIBS = 

#------------------------
# BUILD TOOL COMMANDS
#------------------------

# C Compiler: Lcc-win64 C Compiler
CC_PATH = $(LCC_ROOT)/bin
CC = "$(CC_PATH)/lcc64"

# Linker: Lcc-win64 Linker
LD_PATH = $(LCC_ROOT)/bin
LD = "$(LD_PATH)/lcclnk64"

# Archiver: Lcc-win64 Archiver
AR_PATH = $(LCC_ROOT)/bin
AR = "$(AR_PATH)/lcclib64"

# MEX Tool: MEX Tool
MEX_PATH = $(MATLAB_BIN)
MEX = "$(MEX_PATH)/mex"

# Download: Download
DOWNLOAD =

# Execute: Execute
EXECUTE = $(PRODUCT)

# Builder: GMAKE Utility
MAKE_PATH = %MATLAB%\bin\win64
MAKE = "$(MAKE_PATH)/gmake"


#-------------------------
# Directives/Utilities
#-------------------------

CDEBUG              = -g
C_OUTPUT_FLAG       = -Fo
LDDEBUG             =
OUTPUT_FLAG         = -o
ARDEBUG             =
STATICLIB_OUTPUT_FLAG = /out:
MEX_DEBUG           = -g
RM                  = @del /F
ECHO                = @echo
MV                  = @move
RUN                 =

#----------------------------------------
# "Faster Builds" Build Configuration
#----------------------------------------

ARFLAGS              =
CFLAGS               = -c -w -noregistrylookup -nodeclspec -I$(LCC_ROOT)/include64
DOWNLOAD_FLAGS       =
EXECUTE_FLAGS        =
LDFLAGS              = -s -L$(LCC_LIB) $(LDFLAGS_ADDITIONAL)
MEX_CPPFLAGS         =
MEX_CPPLDFLAGS       =
MEX_CFLAGS           = -win64 $(MEX_SRC) $(MEX_OPT_FILE)$(INCLUDES) -outdir $(RELATIVE_PATH_TO_ANCHOR)
MEX_LDFLAGS          = LINKFLAGS="$$LINKFLAGS $(LDFLAGS_ADDITIONAL)"
MAKE_FLAGS           = -f $(MAKEFILE)
SHAREDLIB_LDFLAGS    = -dll -entry LibMain -s -L$(LCC_LIB) $(LDFLAGS_ADDITIONAL) $(DEF_FILE)

#--------------------
# File extensions
#--------------------

H_EXT               = .h
OBJ_EXT             = .obj
C_EXT               = .c
EXE_EXT             = .exe
SHAREDLIB_EXT       = .dll
STATICLIB_EXT       = .lib
MEX_EXT             = .mexw64
MAKE_EXT            = .mk


###########################################################################
## OUTPUT INFO
###########################################################################

PRODUCT = $(RELATIVE_PATH_TO_ANCHOR)/Calc_pressao.lib
PRODUCT_TYPE = "static-library"
BUILD_TYPE = "Static Library"

###########################################################################
## INCLUDE PATHS
###########################################################################

INCLUDES_BUILDINFO = -I$(START_DIR)/codegen/lib/Calc_pressao -I$(START_DIR) -I$(MATLAB_ROOT)/extern/include -I$(MATLAB_ROOT)/simulink/include -I$(MATLAB_ROOT)/rtw/c/src -I$(MATLAB_ROOT)/rtw/c/src/ext_mode/common -I$(MATLAB_ROOT)/rtw/c/ert

INCLUDES = $(INCLUDES_BUILDINFO)

###########################################################################
## DEFINES
###########################################################################

DEFINES_STANDARD = -DMODEL=Calc_pressao -DHAVESTDIO -DUSE_RTMODEL

DEFINES = $(DEFINES_STANDARD)

###########################################################################
## SOURCE FILES
###########################################################################

SRCS = $(START_DIR)/codegen/lib/Calc_pressao/Calc_pressao_rtwutil.c $(START_DIR)/codegen/lib/Calc_pressao/Calc_pressao_initialize.c $(START_DIR)/codegen/lib/Calc_pressao/Calc_pressao_terminate.c $(START_DIR)/codegen/lib/Calc_pressao/Calc_pressao.c $(START_DIR)/codegen/lib/Calc_pressao/calc_r.c $(START_DIR)/codegen/lib/Calc_pressao/sqrt1.c $(START_DIR)/codegen/lib/Calc_pressao/calc_th.c $(START_DIR)/codegen/lib/Calc_pressao/complexTimes.c $(START_DIR)/codegen/lib/Calc_pressao/calc_v.c $(START_DIR)/codegen/lib/Calc_pressao/calc_H.c $(START_DIR)/codegen/lib/Calc_pressao/fix.c $(START_DIR)/codegen/lib/Calc_pressao/cbesj.c $(START_DIR)/codegen/lib/Calc_pressao/cbinu.c $(START_DIR)/codegen/lib/Calc_pressao/cseri.c $(START_DIR)/codegen/lib/Calc_pressao/log1.c $(START_DIR)/codegen/lib/Calc_pressao/gammaln.c $(START_DIR)/codegen/lib/Calc_pressao/exp.c $(START_DIR)/codegen/lib/Calc_pressao/cuchk.c $(START_DIR)/codegen/lib/Calc_pressao/cmlri.c $(START_DIR)/codegen/lib/Calc_pressao/casyi.c $(START_DIR)/codegen/lib/Calc_pressao/cuoik.c $(START_DIR)/codegen/lib/Calc_pressao/cunhj.c $(START_DIR)/codegen/lib/Calc_pressao/cunik.c $(START_DIR)/codegen/lib/Calc_pressao/cairy.c $(START_DIR)/codegen/lib/Calc_pressao/cbknu.c $(START_DIR)/codegen/lib/Calc_pressao/sinh.c $(START_DIR)/codegen/lib/Calc_pressao/cosh.c $(START_DIR)/codegen/lib/Calc_pressao/calc_p.c $(START_DIR)/codegen/lib/Calc_pressao/sum.c $(START_DIR)/codegen/lib/Calc_pressao/rt_nonfinite.c $(START_DIR)/codegen/lib/Calc_pressao/rtGetNaN.c $(START_DIR)/codegen/lib/Calc_pressao/rtGetInf.c

ALL_SRCS = $(SRCS)

###########################################################################
## OBJECTS
###########################################################################

OBJS = Calc_pressao_rtwutil.obj Calc_pressao_initialize.obj Calc_pressao_terminate.obj Calc_pressao.obj calc_r.obj sqrt1.obj calc_th.obj complexTimes.obj calc_v.obj calc_H.obj fix.obj cbesj.obj cbinu.obj cseri.obj log1.obj gammaln.obj exp.obj cuchk.obj cmlri.obj casyi.obj cuoik.obj cunhj.obj cunik.obj cairy.obj cbknu.obj sinh.obj cosh.obj calc_p.obj sum.obj rt_nonfinite.obj rtGetNaN.obj rtGetInf.obj

ALL_OBJS = $(OBJS)

###########################################################################
## PREBUILT OBJECT FILES
###########################################################################

PREBUILT_OBJS = 

###########################################################################
## LIBRARIES
###########################################################################

LIBS = 

###########################################################################
## SYSTEM LIBRARIES
###########################################################################

SYSTEM_LIBS = 

###########################################################################
## ADDITIONAL TOOLCHAIN FLAGS
###########################################################################

#---------------
# C Compiler
#---------------

CFLAGS_BASIC = $(DEFINES) $(INCLUDES)

CFLAGS += $(CFLAGS_BASIC)

###########################################################################
## INLINED COMMANDS
###########################################################################

###########################################################################
## PHONY TARGETS
###########################################################################

.PHONY : all build clean info prebuild download execute


all : build
	@echo "### Successfully generated all binary outputs."


build : prebuild $(PRODUCT)


prebuild : 


download : build


execute : download


###########################################################################
## FINAL TARGET
###########################################################################

#---------------------------------
# Create a static library         
#---------------------------------

$(PRODUCT) : $(OBJS) $(PREBUILT_OBJS)
	$(PERL) $(GEN_LNK_SCRIPT) $(CMD_FILE) $(subst /,\,$(OBJS))
	@echo "### Creating static library "$(PRODUCT)" ..."
	$(AR) $(ARFLAGS) /out:$(PRODUCT) @$(CMD_FILE)
	@echo "### Created: $(PRODUCT)"
	$(RM) $(CMD_FILE)


###########################################################################
## INTERMEDIATE TARGETS
###########################################################################

#---------------------
# SOURCE-TO-OBJECT
#---------------------

%.obj : %.c
	$(CC) $(CFLAGS) -Fo"$@" $(subst /,\,"$<")


%.obj : $(RELATIVE_PATH_TO_ANCHOR)/%.c
	$(CC) $(CFLAGS) -Fo"$@" $(subst /,\,"$<")


%.obj : $(START_DIR)/%.c
	$(CC) $(CFLAGS) -Fo"$@" $(subst /,\,"$<")


%.obj : $(START_DIR)/codegen/lib/Calc_pressao/%.c
	$(CC) $(CFLAGS) -Fo"$@" $(subst /,\,"$<")


%.obj : $(MATLAB_ROOT)/rtw/c/src/%.c
	$(CC) $(CFLAGS) -Fo"$@" $(subst /,\,"$<")


###########################################################################
## DEPENDENCIES
###########################################################################

$(ALL_OBJS) : $(MAKEFILE) rtw_proj.tmw


###########################################################################
## MISCELLANEOUS TARGETS
###########################################################################

info : 
	@echo "### PRODUCT = $(PRODUCT)"
	@echo "### PRODUCT_TYPE = $(PRODUCT_TYPE)"
	@echo "### BUILD_TYPE = $(BUILD_TYPE)"
	@echo "### INCLUDES = $(INCLUDES)"
	@echo "### DEFINES = $(DEFINES)"
	@echo "### ALL_SRCS = $(ALL_SRCS)"
	@echo "### ALL_OBJS = $(ALL_OBJS)"
	@echo "### LIBS = $(LIBS)"
	@echo "### MODELREF_LIBS = $(MODELREF_LIBS)"
	@echo "### SYSTEM_LIBS = $(SYSTEM_LIBS)"
	@echo "### TOOLCHAIN_LIBS = $(TOOLCHAIN_LIBS)"
	@echo "### CFLAGS = $(CFLAGS)"
	@echo "### LDFLAGS = $(LDFLAGS)"
	@echo "### SHAREDLIB_LDFLAGS = $(SHAREDLIB_LDFLAGS)"
	@echo "### ARFLAGS = $(ARFLAGS)"
	@echo "### MEX_CFLAGS = $(MEX_CFLAGS)"
	@echo "### MEX_CPPFLAGS = $(MEX_CPPFLAGS)"
	@echo "### MEX_LDFLAGS = $(MEX_LDFLAGS)"
	@echo "### MEX_CPPLDFLAGS = $(MEX_CPPLDFLAGS)"
	@echo "### DOWNLOAD_FLAGS = $(DOWNLOAD_FLAGS)"
	@echo "### EXECUTE_FLAGS = $(EXECUTE_FLAGS)"
	@echo "### MAKE_FLAGS = $(MAKE_FLAGS)"


clean : 
	$(ECHO) "### Deleting all derived files..."
	$(RM) $(subst /,\,$(PRODUCT))
	$(RM) $(subst /,\,$(ALL_OBJS))
	$(ECHO) "### Deleted all derived files."


