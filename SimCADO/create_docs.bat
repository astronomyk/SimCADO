cd ./simcado/docs/
python Keywords.py
sphinx-apidoc --force -M --separate --no-toc -o source/API ../
sphinx-build -b html -a -j 6 ./source ./build
pause