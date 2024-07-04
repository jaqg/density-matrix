# #!/bin/sh
# +--------------------------------------------+
# | Author: Jose Antonio Quinonero Gris        |
# | Creation date: Tuesday 19:12:57 25/06/2024 |
# +--------------------------------------------+

# Verify that three arguments are provided
# if [ "$#" -ne 3 ]; then
#     echo "Usage: $0 molecule n level"
#     exit 1
# fi

# Get the arguments
# molecule="$1"
# n="$2"
# level="$3"
molecule="H2O"
n=8
level="MCSCF"

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

# Initialize all_energies and all_errors files if they do not exist
if [ ! -f "$DIR/all_energies" ]; then
    echo -n "Functional" > "$DIR/all_energies"
fi

if [ ! -f "$DIR/all_errors" ]; then
    echo -n "Functional" > "$DIR/all_errors"
fi

# Run ./main for each molecule-i-level.mol file
for i in $(seq 1 "$n"); do
    OUTPUT="$DIR/INPUT"
    while IFS= read -r line; do
        case "$line" in
            "folder") echo "data" ;;
            "SIRIFC") echo "${molecule}-${i}-${level}-SIRIFC" ;;
            "molecule.fmt") echo "${molecule}-${i}-${level}-integrals" ;;
            *) echo "$line" ;;
        esac
    done < "$TEMPLATE" > "$OUTPUT"

    # Run the program ./main
    ./main

    # Process the ENERGIES file and update all_energies and all_errors
    if [ -f "ENERGIES" ]; then
        if [ "$i" -eq 1 ]; then
            # For the first iteration, add the headers to all_energies and all_errors
            paste "$DIR/all_energies" <(echo -e "\n${molecule}-${i}-${level}") > "$DIR/all_energies.tmp" && mv "$DIR/all_energies.tmp" "$DIR/all_energies"
            paste "$DIR/all_errors" <(echo -e "\n${molecule}-${i}-${level}") > "$DIR/all_errors.tmp" && mv "$DIR/all_errors.tmp" "$DIR/all_errors"
        else
            # For subsequent iterations, add the i value to the first row
            sed -i "1s/$/\t${i}/" "$DIR/all_energies"
            sed -i "1s/$/\t${i}/" "$DIR/all_errors"
        fi
        # Add the second column of ENERGIES to all_energies
        paste "$DIR/all_energies" <(cut -d' ' -f2 ENERGIES) > "$DIR/all_energies.tmp" && mv "$DIR/all_energies.tmp" "$DIR/all_energies"
        # Add the third column of ENERGIES to all_errors
        paste "$DIR/all_errors" <(cut -d' ' -f3 ENERGIES) > "$DIR/all_errors.tmp" && mv "$DIR/all_errors.tmp" "$DIR/all_errors"
    else
        echo "Error: the ENERGIES file is not found after running ./main."
        exit 1
    fi

    # Success message
    echo "Successfull execution for ${molecule}-${i}-${level}."
done

