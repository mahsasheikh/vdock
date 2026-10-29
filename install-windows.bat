@echo off
REM Double-click or run from a terminal: installs the vdock environment on Windows.
cd /d "%~dp0"
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0install-windows.ps1" %*
pause
