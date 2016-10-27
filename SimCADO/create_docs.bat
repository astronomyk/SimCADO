sphinx-apidoc --full --force -o simcado/docs/source/API simcado
cd simcado/docs/
sphinx-build -b html source/API site/API/_build
mkdocs build --clean
cd ../..
