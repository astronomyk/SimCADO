D:
cd D:\Dropbox\Uni\PhD\SimCADO\SimCADO

REM #######################################
REM Compile the new distribution
REM #######################################

del .\dist\*.tar.gz
del .\dist\*.zip
dir .\dist\
python setup.py sdist

@echo off

setlocal enabledelayedexpansion
set params=
for /f "delims=" %%a in ('dir .\dist\* /s/b') do set params=!params! %%a 
echo !params!

pip install -I !params!

REM #######################################
REM Run tests on the installed distribution
REM #######################################

cd ../dist_tests
REM pytest test_run.py
python test_play.py
cd ../SimCADO

REM #######################################
REM Make the f*ing docs
REM #######################################

create_docs.bat

echo apple is not the hammer 
pause
