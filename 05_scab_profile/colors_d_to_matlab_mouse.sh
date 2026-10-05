#!/bin/bash

input_file="$1"
output_file="${input_file}.txt"

sed 's/chr//g; s/^X	/20	/g; s/^Y	/21	/g' "${input_file}" > "${output_file}"

