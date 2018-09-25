cd src

find energy/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format-3.5 -i  -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0}"
find histogram/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format-3.5 -i -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0}"
find include/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format-3.5 -i -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0}"
find io/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format-3.5 -i -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0}"
find main/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format-3.5 -i -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0}"
find mc/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format-3.5 -i -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0}"
find mersenne/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format-3.5 -i -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0}"
find polarization/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format-3.5 -i -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0}"
find quantum_rotation/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format-3.5 -i -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0}"

cd ..
