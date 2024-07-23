#!/bin/bash

# Function to extract the total energy from a DALTON output file
extract_energy() {
    local output_file=$1
    local energy_keyword="Final MCSCF energy"

    # Extract the line containing the energy keyword
    energy_line=$(grep "$energy_keyword" "$output_file" | tail -n 1)

    # Extract the energy value
    energy_value=$(echo "$energy_line" | awk -F: '{print $2}')

    echo "$energy_value"
}

# File to store all energies
summary_file="reference_energies.dat"
echo "File Energy" > "$summary_file"

# Loop through each *.out file in the current directory
for output_file in *.out; do
    # Extract the base filename without the .out extension
    system=$(basename "$output_file" .out)

    # Extract the last energy from the output file
    energy=$(extract_energy "$output_file")
    
    # Append the filename and energy to the summary file
    echo "$system $energy" >> "$summary_file"
done
