#!/bin/bash

# Directory containing the header files
directory="."

# Line to look for and the line to insert
search_line="#define ALLOW_SURFVEL_UNCERT"
insert_line="#define SURVEL_UNCERT_FIELD 'vs_uncert'"

# Loop through all files in the directory
for file in "$directory"/*.h; do
  if grep -q "$search_line" "$file"; then
    # Use sed to insert the line after the search line
    sed -i "/$search_line/a \\$insert_line" "$file"
  fi
done

echo "Header files updated successfully."
