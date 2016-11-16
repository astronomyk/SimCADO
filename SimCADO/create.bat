D:
cd D:\Dropbox\Uni\PhD\SimCADO\SimCADO

del .\dist\*.zip
dir .\dist\
python setup.py sdist

@echo off

setlocal enabledelayedexpansion
set params=
for /f "delims=" %%a in ('dir .\dist\*.zip /s/b') do set params=!params! %%a 
echo !params!

pip install -I !params!

echo apple is not the hammer 
pause
