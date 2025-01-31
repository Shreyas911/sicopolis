#!/bin/bash

# Define the directory containing Fortran header files
DIR="."

# Define the search pattern (exact match)
SEARCH="#define NUM_CTRL_GENARR2D 13"

# Define the replacement text with properly escaped backslashes
REPLACE='#define NUM_CTRL_GENARR2D 17'

# Process all .h files in the directory
for file in "$DIR"/*.h; do
    if [[ -f "$file" ]]; then
        # Use sed to replace the text and create a backup
        sed -i.bak -e "/$SEARCH/c\\
$REPLACE" "$file"
        echo "Processed: $file"
    fi
done

echo "Replacement completed."

