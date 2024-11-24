#!/bin/bash

# Default values
n=10
PROGRAM_PATH="./my_program"

# Function to display usage
usage() {
    echo "Usage: $0 [-n number_of_runs] [-p program_path]"
    exit 1
}

# Parse command-line arguments
while getopts ":n:p:" opt; do
  case $opt in
    n)
      n=$OPTARG
      ;;
    p)
      PROGRAM_PATH=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      ;;
  esac
done

# Check if the program exists and is executable
#if ! [ -x "$PROGRAM_PATH" ]; then
#    echo "Error: Program '$PROGRAM_PATH' not found or is not executable."
#    exit 1
#fi

# Loop to run the program 'n' times
for ((i=1; i<=n; i++))
do
    # Run the program and capture its output
    output=$($PROGRAM_PATH)

    # Extract the 'Time taken' value using awk
    time_taken=$(echo "$output" | awk -F'Time taken: ' '{print $2}' | awk '{print $1}' | sed 's/ms//')

    # Display the extracted 'Time taken' value
    echo "${time_taken}"
done
