python Keywords.py
sphinx-apidoc -f -M -e -T -o source/API ../
sphinx-build -b html -a -j 6 ./source ./build
pause