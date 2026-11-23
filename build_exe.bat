@echo off
cd /d "%~dp0"
echo ============================================
echo   vanVleck - build Windows executable
echo ============================================
echo [1/3] Installing dependencies (PySide6, numpy, PyInstaller)...
python -m pip install --upgrade -r requirements.txt pyinstaller
if errorlevel 1 goto err
echo [2/3] Building vanVleck.exe (this takes 1-3 minutes)...
python -m PyInstaller --noconfirm --clean --onefile --windowed --name vanVleck main.py
if errorlevel 1 goto err
copy /y dist\vanVleck.exe vanVleck.exe >nul
echo [3/3] Done!  %~dp0vanVleck.exe
echo You can double-click vanVleck.exe, or drag a .cif/.xyz file onto it.
pause
exit /b 0
:err
echo.
echo Build FAILED. Make sure Python is installed and "python" works in cmd.
pause
exit /b 1
