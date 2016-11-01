D:
cd D:\Dropbox\Uni\PhD\SimCADO\SimCADO

cd simcado/docs/
mkdocs build --clean

cd site/API/
mkdir _source
mkdir _build

cd ../../../..
sphinx-apidoc -f -F -o simcado/docs/site/API/_source simcado

cd simcado/docs/site/API/_source
copy ..\..\..\conf.py conf.py
sphinx-build -b html -E .\ ..\_build