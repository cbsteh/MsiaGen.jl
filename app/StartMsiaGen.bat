@echo off
rem Double-click to start MsiaGen on Windows. Keep this window open while using it.
cd /d "%~dp0"
cls
echo.
echo   MSIAGEN is starting. Please wait, this takes about 1 minute.
echo   Your web browser will open by itself when it is ready.
echo.
julia --project=. -t auto start.jl
pause
