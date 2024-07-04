#!/bin/sh

dir='.'

find "$dir" -type f -name "*.tex" -print0 | while IFS= read -r -d '' file; do
    # cat "$file" | aspell list | sort -u
    echo "$file" ; aspell list < "$file" | sort | uniq -c
done
