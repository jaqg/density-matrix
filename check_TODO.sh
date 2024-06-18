#!/bin/sh

# grep recursively for TODO in the directory 'directory'
directory='.'
script_name=$(basename "$0")

# skip direcctories that start with '.' and this script
find "$directory" -type d -name '.git' -prune -o -type f ! -name "$script_name" -exec grep -Hn "TODO" {} \;
