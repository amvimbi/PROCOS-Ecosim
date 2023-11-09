@echo off
rem Make sure VS scripts find the reg command:
set PATH=C:\Windows\System32;C:\PROGRA~1\ECOSIM~1\ECOSIM~1.2\platform\WIN64_~1;
rem echo "%1"
if "%1" == "" (goto skip_detection)
rem Call compiler script if possible ...
if "%1" == "win64_vc2017" (goto find_vc2017)
if "%1" == "win64_vc2017icc14" (goto find_vc2017)
:find_vc2017
set COMPILER_DIR_VC2017=
for /f "usebackq tokens=1* delims=: " %%i in (`C:\PROGRA~1\ECOSIM~1\ECOSIM~1.2\platform\WIN64_~1\vswhere.exe -latest -requires Microsoft.VisualStudio.Workload.NativeDesktop`) do (if /i "%%i" == "installationPath" set COMPILER_DIR_VC2017=%%j)
:continue_configuring_compiler
set COMPILER_SCRIPT=
if "%1" == "win32_vc2008" (set COMPILER_SCRIPT="%VS90COMNTOOLS%vsvars32.bat")
if "%1" == "win32_vc2010" (set COMPILER_SCRIPT="%VS100COMNTOOLS%vsvars32.bat")
if "%1" == "win64_vc2010" (set COMPILER_SCRIPT="%VS100COMNTOOLS%\..\..\VC\bin\amd64\vcvars64.bat")
if "%1" == "win64_vc2010icc14" (set COMPILER_SCRIPT="%ICPP_COMPILER14%bin\iclvars.bat")
if "%1" == "win64_vc2013" (set COMPILER_SCRIPT="%VS120COMNTOOLS%\..\..\VC\bin\amd64\vcvars64.bat")
if "%1" == "win64_vc2013icc14" (set COMPILER_SCRIPT="%ICPP_COMPILER14%bin\iclvars.bat")
if "%1" == "win64_vc2015" (set COMPILER_SCRIPT="%VS140COMNTOOLS%\..\..\VC\bin\amd64\vcvars64.bat")
if "%1" == "win64_vc2015icc14" (set COMPILER_SCRIPT="%ICPP_COMPILER14%bin\iclvars.bat")
if "%1" == "win64_vc2017" (set COMPILER_SCRIPT="%COMPILER_DIR_VC2017%\VC\Auxiliary\Build\vcvars64.bat")
if "%1" == "win64_vc2017icc14" (set COMPILER_SCRIPT="%ICPP_COMPILER14%bin\iclvars.bat")
echo Calling compiler script (if necessary)...
if "%1" == "win32_vc2008" (call %COMPILER_SCRIPT%)
if "%1" == "win32_vc2010" (call %COMPILER_SCRIPT%)
if "%1" == "win64_vc2010" (call %COMPILER_SCRIPT%)
if "%1" == "win64_vc2010icc14" (call %COMPILER_SCRIPT% intel64 vs2010)
if "%1" == "win64_vc2013" (call %COMPILER_SCRIPT%)
if "%1" == "win64_vc2013icc14" (call %COMPILER_SCRIPT% intel64 vs2013)
if "%1" == "win64_vc2015" (call %COMPILER_SCRIPT%)
if "%1" == "win64_vc2015icc14" (call %COMPILER_SCRIPT% intel64 vs2013)
if "%1" == "win64_vc2017" (call %COMPILER_SCRIPT%)
if "%1" == "win64_vc2017icc14" (call %COMPILER_SCRIPT% intel64 vs2013)
if "%1" == "win32_gcc32v4_4" goto set_win32_gcc32v4_4
if "%1" == "win32_gcc32v4_7" goto set_win32_gcc32v4_7
if "%1" == "win64_gcc64v4_9" goto set_win64_gcc64v4_9
goto quit
:set_win32_gcc32v4_4
set PATH=C:\Program Files\EcosimPro\EcosimPro_6.0.2\3rdparty\compiler\mingw_gcc_4_4_0\MinGW\bin
set INCLUDE=
set LIB=
goto quit
:set_win32_gcc32v4_7
set PATH=C:\Program Files\EcosimPro\EcosimPro_6.0.2\3rdparty\compiler\mingw_gcc_4_7\MinGW\bin
set INCLUDE=
set LIB=
goto quit
:set_win64_gcc64v4_9
set PATH=C:\Program Files\EcosimPro\EcosimPro_6.0.2\3rdparty\compiler\mingw64_gcc_4_9\mingw64\bin
set INCLUDE=
set LIB=
goto quit
:skip_detection
:quit
nmake all /f "C:/Program Files/EcosimPro/EcosimPro_6.0.2/mk/makefile_default_win64_vc.mak"