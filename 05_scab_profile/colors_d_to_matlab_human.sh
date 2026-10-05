#!/bin/bash

input_file="$1"
output_file="${input_file}.txt"

sed 's/^X	/23	/g; s/^Y	/24	/g' "${input_file}" > "${output_file}"
