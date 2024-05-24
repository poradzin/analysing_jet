#!/bin/bash

# Check if the correct number of arguments is provided
if [[ $# -lt 5 || $# -gt 6 ]]; then
    echo "Usage: $0 <directory> <section> <letter> <start_number> <end_number> [w]"
    exit 1
fi

# Read the arguments
directory=$1
section=$2
letter=$3
start_number=$(printf "%02d" $4)
end_number=$(printf "%02d" $5)

# Check if the script should write to a file
write_to_file=false
if [[ $# -eq 6 && $6 == "w" ]]; then
    write_to_file=true
fi

# Specify the starting directory
start_dir="/common/transp_shared/Data/result/JET/$directory"

# Create the output file name based on the provided section and range
output_file="${directory}_${section}_${letter}${start_number}_to_${letter}${end_number}.txt"

# Clear the output file if it exists and writing is enabled
if [[ $write_to_file == true ]]; then
    > $output_file
fi

# Loop through the specified range of directories
for i in $(seq -f "%02g" $4 $5); do
    dir="$start_dir/${letter}$i"
    file="$dir/${directory}${letter}${i}_Provenance.DAT"
    
    # Check if the file exists
    if [[ -f "$file" ]]; then
        # Extract lines from the specified section until a blank line is encountered
        sec_lines=$(awk -v sec="$section" '$0 ~ "^"sec":" {flag=1; next} flag && NF {print; next} flag && !NF {flag=0}' "$file")

        # Extract DDA, UID, and SEQ values from the specified section
        DDA=$(echo "$sec_lines" | grep "^DDA" | awk -F ' : ' '{print $2}')
        UID_VAL=$(echo "$sec_lines" | grep "^UID" | awk -F ' : ' '{print $2}')
        SEQ=$(echo "$sec_lines" | grep "^SEQ" | awk -F ' : ' '{print $2}')

        # Print the result to the terminal
        echo "${letter}$i:$section: $DDA/$UID_VAL/$SEQ"
        
        # Save the result to the output file if writing is enabled
        if [[ $write_to_file == true ]]; then
            echo "${letter}$i:$section: $DDA/$UID_VAL/$SEQ" >> $output_file
        fi
    else
        echo "Run ${letter}${i} does not exist"  # Optional: handle missing files
    fi
done

