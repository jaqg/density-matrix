# #!/bin/sh
# +--------------------------------------------+
# | Author: Jose Antonio Quinonero Gris        |
# | Creation date: Tuesday 19:12:57 25/06/2024 |
# +--------------------------------------------+

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 data_folder level"
    exit 1
fi

# Get the arguments
# molecule="$1"
# n="$2"
molecule="H2O"
n=8
# level="MCSCF"
# data_folder="data"
data_folder="$1"
level="$2"

# Verify that n is a positive integer
if ! [[ "$n" =~ ^[0-9]+$ ]]; then
    echo "Error: n must be a positive integer."
    exit 1
fi

# Get the directory of the script
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Full path to the input_template file
TEMPLATE="$DIR/input_template"

# Check if the input_template file exists
if [ ! -f "$TEMPLATE" ]; then
    echo "Error: the input_template file does not exist in the script directory."
    exit 1
fi

# Initialize all_energies, all_errors, all_minkowski
echo -e "#	LS	       BBC1	      BBC2	     BBC3	    BBC3M" > all_energies
echo -e "#	LS	     BBC1	  BBC2	   BBC3	    BBC3M" > all_errors
echo -e "#	LS	     BBC1	  BBC2	   BBC3	    BBC3M" > all_minkowski

# Run ./main for each molecule-i-level.mol file
for i in $(seq 1 "$n"); do
    OUTPUT="$DIR/INPUT"
    while IFS= read -r line; do
        case "$line" in
            "folder") echo "$data_folder" ;;
            "SIRIFC") echo "${molecule}-${i}-${level}-SIRIFC" ;;
            "molecule.fmt") echo "${molecule}-${i}-${level}-integrals" ;;
            *) echo "$line" ;;
        esac
    done < "$TEMPLATE" > "$OUTPUT"

    # Run the program ./main
    ./main

    # Process the plot_data file -> all_energies, all_errors, all_minkowski
    if [ -f "plot_data" ]; then
        
        # Extract energies, errors and minkowski distances from plot_data
        energies=$(grep "E (Ha) :" plot_data | awk '{print $4, $5, $6, $7, $8}')
        errors=$(grep "ΔE (Ha):" plot_data | awk '{print $3, $4, $5, $6, $7}' | sed 's/D/e/g')
        minkowski=$(grep "Mink. d:" plot_data | awk '{print $3, $4, $5, $6, $7}')

        # Append the extracted data to the respective files with the index i
        echo -e "$i\t$energies" >> all_energies
        echo -e "$i\t$errors" >> all_errors
        echo -e "$i\t$minkowski" >> all_minkowski
    fi

    # Success message
    echo "Successfull execution for ${molecule}-${i}-${level}."
done

cp all_energies  "$data_folder"
cp all_errors    "$data_folder"
cp all_minkowski "$data_folder"
