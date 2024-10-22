@echo off
set MATLAB=C:\PROGRA~1\MATLAB\R2017a
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\R2017a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=Calc_pressao_2_mex
set MEX_NAME=Calc_pressao_2_mex
set MEX_EXT=.mexw64
call "C:\PROGRA~1\MATLAB\R2017a\sys\lcc64\lcc64\mex\lcc64opts.bat"
echo # Make settings for Calc_pressao > Calc_pressao_mex.mki
echo COMPILER=%COMPILER%>> Calc_pressao_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> Calc_pressao_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> Calc_pressao_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> Calc_pressao_mex.mki
echo LINKER=%LINKER%>> Calc_pressao_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> Calc_pressao_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> Calc_pressao_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> Calc_pressao_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> Calc_pressao_mex.mki
echo BORLAND=%BORLAND%>> Calc_pressao_mex.mki
echo OMPFLAGS= >> Calc_pressao_mex.mki
echo OMPLINKFLAGS= >> Calc_pressao_mex.mki
echo EMC_COMPILER=lcc64>> Calc_pressao_mex.mki
echo EMC_CONFIG=optim>> Calc_pressao_mex.mki
"C:\Program Files\MATLAB\R2017a\bin\win64\gmake" -B -f Calc_pressao_mex.mk
