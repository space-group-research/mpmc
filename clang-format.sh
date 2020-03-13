cd src

find energy/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format -i  -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0, SortIncludes: false}"
find histogram/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0, SortIncludes: false}"
find include/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0, SortIncludes: false}"
find io/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0, SortIncludes: false}"
find main/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0, SortIncludes: false}"
find mc/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0, SortIncludes: false}"
find polarization/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0, SortIncludes: false}"
find quantum_rotation/ -iname *.h -o -iname *.c -o -iname *.cpp | xargs clang-format -i -style="{BasedOnStyle: google, IndentWidth: 4, ColumnLimit: 0, SortIncludes: false}"

cd ..

