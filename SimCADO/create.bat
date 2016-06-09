del .\dist\*.zip
dir .\dist\
python setup.py sdist
pip uninstall simcado
pip install -I .\dist\simcado-0.2dev.zip

@echo off
echo apple is not the hammer 
pause
