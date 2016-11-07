D:
cd D:\Dropbox\Uni\PhD\SimCADO\SimCADO

cd simcado/docs/source/
mkdir images
cd ../

python SimCADO_defaults.py
mkdocs build --clean

cd site/SimCADO_defaults/
mkdir images
copy ..\..\source\images\*.png .\images\

cd ../API/
mkdir _source
mkdir _build


cd ../../../..
sphinx-apidoc -f -F -o simcado/docs/site/API/_source simcado

cd simcado/docs/site/API/_source
copy ..\..\..\conf.py conf.py
sphinx-build -b html -E .\ ..\_build

pause