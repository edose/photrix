echo off
set "START_UTC=2"
set "STOP_UTC=12"
echo 
C:\Programs\Python_envs\py_3.7.5_venv\Scripts\python.exe C:/Dev/photrix/photrix/fits_monitor.py %1 %START_UTC% %STOP_UTC%
pause