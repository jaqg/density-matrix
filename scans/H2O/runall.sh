#!/bin/sh
# +-------------------------------------------+
# | Author: Jose Antonio Quinonero Gris       |
# | Creation date: Sunday 23:18:50 16/06/2024 |
# +-------------------------------------------+

# Original input file (template)
input_template="input_template"

# Number of files
num_files=8

# Read the template input file
template=$(cat $input_template)

# Initialize the all_energies.dat file with the header
echo -n "#  " > all_energies.dat
head -n 1 energies.dat | awk '{print $0}' >> all_energies.dat

# Rewrite 'input' for each set of H2O-n-* files
for n in $(seq 0 $((num_files - 1)))
do
    # Replace H2O-0- with H2O-n-
    input_content=$(echo "$template" | sed "s/H2O-0-/H2O-$n-/g")
    
    # Save the updated content to the input file
    echo "$input_content" > input
    
    # Run the 'main' script
    ./main

    # Check if 'energies.dat' exists and has the expected format
    if [ -f energies.dat ]; then
        # Read the second line and append it to 'all_energies.dat'
        echo -n "$n " >> all_energies.dat
        sed -n '2p' energies.dat >> all_energies.dat
    else
        echo "Warning: energies.dat not found or empty for iteration $n"
    fi
    
    echo "Executed main with input for H2O-$n"
done
