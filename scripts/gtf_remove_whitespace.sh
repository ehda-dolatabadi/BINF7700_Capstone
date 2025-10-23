#!/usr/bin/env bash
set -euo pipefail

in="$1"
out="${in%.gtf}.noWS.gtf"

awk '{
  s=$0; out="";
  while (match(s, /"[^"]*"/)) {
    out = out substr(s,1,RSTART-1);
    q = substr(s,RSTART+1,RLENGTH-2);
    gsub(/[ \t\r\f\v]+/, "", q);      # remove spaces/tabs/etc inside quotes
    out = out "\"" q "\"";
    s = substr(s,RSTART+RLENGTH);
  }
  print out s
}' "$in" > "$out"

echo "wrote: $out"
