#!/usr/bin/env bash
# Script: gtf_remove_whitespace.sh
# Purpose: Remove whitespace from GTF file attributes while preserving structure
# Usage: ./gtf_remove_whitespace.sh <input.gtf>

# Exit on error, undefined variables, and pipe failures
set -euo pipefail

# Input GTF file from command line argument
in="$1"

# Output file: original filename with .noWS.gtf extension
out="${in%.gtf}.noWS.gtf"

# Use AWK to process GTF file line by line
# Strategy: Find quoted strings, remove whitespace within them, preserve structure
 awk '{
  s=$0; out="";
  # Loop through all quoted strings in the line
  while (match(s, /"[^"]*"/)) {
    # Keep everything before the quote
    out = out substr(s,1,RSTART-1);
    # Extract the quoted content (without quotes)
    q = substr(s,RSTART+1,RLENGTH-2);
    # Remove all whitespace characters (spaces, tabs, etc.) from quoted content
    gsub(/[ \t\r\f\v]+/, "", q);
    # Add back the quotes around the cleaned content
    out = out "\"" q "\"";
    # Continue with the rest of the line
    s = substr(s,RSTART+RLENGTH);
  }
  # Print the processed line
  print out s
}' "$in" > "$out"

echo "wrote: $out"
