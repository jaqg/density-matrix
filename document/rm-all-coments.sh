#!/bin/bash

# Function to process a .tex file
process_file() {
    local file="$1"
    local tempfile=$(mktemp)  # Create a temporary file

    # Check if the file is newcommands.tex; if yes, skip processing
    if [[ "$file" == *"newcommands.tex" ]]; then
        echo "Skipping $file"
        return
    fi

    # Pattern to match '\begin{comentario}' or '\end{comentario}'
    pattern='\\(begin|end){comentario}'

    # Process each line of the file
    while IFS= read -r line; do
        if [[ "$line" =~ $pattern ]]; then
            # Skip lines containing '\begin{comentario}' or '\end{comentario}'
            continue
        fi
        echo "$line" >> "$tempfile"
    done < "$file"

    # Move the temporary file back to the original file
    mv "$tempfile" "$file"
}

# Function to search for .tex files recursively and process them
search_and_process() {
    local dir="$1"

    # Find .tex files recursively
    find "$dir" -type f -name "*.tex" -print0 | while IFS= read -r -d '' file; do
        process_file "$file"
    done
}

# Check arguments
if [ $# -eq 0 ]; then
    echo "Usage: $0 directory"
    exit 1
fi

# Base directory
base_dir="$1"

# Check if the directory exists
if [ ! -d "$base_dir" ]; then
    echo "Directory '$base_dir' does not exist."
    exit 1
fi

# Call the function to search and process .tex files
search_and_process "$base_dir"

echo "Processing completed."
