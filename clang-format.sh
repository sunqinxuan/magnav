find ./src     -name *.cpp | xargs clang-format -i -style=file
find ./src -name *.hpp | xargs clang-format -i -style=file
find ./src -name *.h | xargs clang-format -i -style=file
