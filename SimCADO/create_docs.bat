cd simcado/docs/
mkdocs build --clean

cd site/API/
mkdir _source
mkdir _build

cd _source
sphinx-apidoc -f -F -o .\ simcado
copy ..\..\..\conf.py conf.py

sphinx-build -b html .\ ..\_build

cd ..\..\..\..\..
pause