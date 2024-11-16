#!/bin/sh  
set -e 
for file in *.bed; do  
    awk '{print "chr" $0}' "$file" > "${file}.tmp" && mv "${file}.tmp" "$file"  
done
echo "done"

