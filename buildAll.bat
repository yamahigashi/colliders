@echo off
setlocal

set "COMPILER=Visual Studio 17 2022"
set "SUPPORTED_MAYA_VERSIONS=2022 2023 2024 2025 2026 2027"
for %%I in ("%~dp0.") do set "SOURCE_DIR=%%~fI"

if not "%~2"=="" goto :usage
if "%~1"=="" goto :buildAll

for %%V in (%SUPPORTED_MAYA_VERSIONS%) do (
    if "%~1"=="%%V" goto :buildOne
)

echo ERROR: Unsupported Maya version "%~1".
goto :usage

:buildOne
call :buildForMayaVersion %~1
if errorlevel 1 (
    echo.
    echo ERROR: Failed to build yddColliders for Maya %~1.
    exit /b 1
)

echo.
echo Successfully built and deployed Maya %~1 to "%SOURCE_DIR%\plugins\%~1".
exit /b 0

:buildAll
for %%V in (%SUPPORTED_MAYA_VERSIONS%) do (
    call :buildForMayaVersion %%V
    if errorlevel 1 (
        echo.
        echo ERROR: Failed to build yddColliders for Maya %%V.
        exit /b 1
    )
)

echo.
echo Successfully built and deployed all Maya versions to "%SOURCE_DIR%\plugins".
exit /b 0

:usage
echo Usage: %~nx0 [Maya version]
echo Supported Maya versions: %SUPPORTED_MAYA_VERSIONS%
exit /b 2

:buildForMayaVersion
setlocal
set "MAYA_VERSION=%~1"
set "BUILD_DIR=%SOURCE_DIR%\mayabuild_%MAYA_VERSION%"
set "PLUGIN_DIR=%SOURCE_DIR%\plugins\%MAYA_VERSION%"
set "MAYA_INSTALL_DIR="

echo.
echo ============================================================
echo Building yddColliders for Maya %MAYA_VERSION%
echo ============================================================

for /f "skip=2 tokens=2*" %%A in ('reg query "HKEY_LOCAL_MACHINE\SOFTWARE\Autodesk\Maya\%MAYA_VERSION%\Setup\InstallPath" /v "MAYA_INSTALL_LOCATION" 2^>nul') do set "MAYA_INSTALL_DIR=%%B"

if not defined MAYA_INSTALL_DIR (
    echo ERROR: Maya %MAYA_VERSION% is not installed or its registry entry is missing.
    exit /b 1
)

for %%I in ("%MAYA_INSTALL_DIR%.") do set "MAYA_DEVKIT_DIR=%%~fI"
if not exist "%MAYA_DEVKIT_DIR%\include" (
    echo ERROR: Maya headers were not found at "%MAYA_DEVKIT_DIR%\include".
    exit /b 1
)
if not exist "%MAYA_DEVKIT_DIR%\lib" (
    echo ERROR: Maya libraries were not found at "%MAYA_DEVKIT_DIR%\lib".
    exit /b 1
)

cmake ^
    -S "%SOURCE_DIR%" ^
    -B "%BUILD_DIR%" ^
    -G "%COMPILER%" ^
    -DMAYA_VERSION=%MAYA_VERSION% ^
    -DMAYA_DEVKIT_DIR="%MAYA_DEVKIT_DIR%"
if errorlevel 1 exit /b 1

cmake --build "%BUILD_DIR%" --config Release
if errorlevel 1 exit /b 1

if not exist "%BUILD_DIR%\Release\yddColliders.mll" (
    echo ERROR: Build succeeded but yddColliders.mll was not produced.
    exit /b 1
)

if not exist "%PLUGIN_DIR%" mkdir "%PLUGIN_DIR%"
copy /y "%BUILD_DIR%\Release\yddColliders.mll" "%PLUGIN_DIR%\yddColliders.mll" >nul
if errorlevel 1 exit /b 1

echo Deployed "%PLUGIN_DIR%\yddColliders.mll".
exit /b 0
