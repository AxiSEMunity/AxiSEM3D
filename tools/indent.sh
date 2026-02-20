#!/bin/sh
#
# This script indents the source code using clang-format.
#
# Usage:
#   ./tools/indent.sh

set -e

# check we are in the project root directory:
if [ ! -f "CMakeLists.txt" ]; then
    echo "Error: not in the project root directory."
    echo "Please run this script from the project root directory."
    exit 1
fi

# check if clang-format is installed:
if ! command -v clang-format > /dev/null 2>&1; then
    echo "Error: clang-format could not be found."
    echo "Please install clang-format 21.1, for example by using the conda environment."
    exit 1
fi

# find all source code files:
source_files=$(find ./src/ -name "*.cpp" -o -name "*.hpp" -o -name "*.h" -o -name "*.c")

# indent the source code:
for file in $source_files; do
    tmpname=$file.clangformat.tmp
    clang-format $file > $tmpname 2>&1

    # overwrite the file with the indented version if the file is different from the original:
    if diff $file $tmpname > /dev/null 2>&1; then
        rm "$tmpname"
    else
        echo "Indented $file"
        mv "$tmpname" "$file"
    fi
done
