#!/usr/bin/env bash
# Usage: bash extract_contacts_unisex.sh /path/to/processed /path/to/folders.txt
# Optional third argument: output.tsv
set -euo pipefail

processed=${1:-./processed}
folder_list=${2:-./folders.txt}
output=${3:-$HOME/research/plateC/contacts_unisex_summary.tsv}

[[ -d "$processed" ]] || { printf 'Directory not found: %s\n' "$processed" >&2; exit 1; }
[[ -r "$folder_list" ]] || { printf 'Folder list not readable: %s\n' "$folder_list" >&2; exit 1; }
mkdir -p -- "$(dirname -- "$output")"
tmp=$(mktemp "${output}.tmp.XXXXXX")
issues=$(mktemp "${output}.issues.tmp.XXXXXX")
trap 'rm -f -- "$tmp" "$issues"' EXIT
printf 'folder_name\tsample_id\tvalue_1\tvalue_2\tvalue_3\tvalue_4\tR1_read_count\n' > "$tmp"
printf 'folder_name\tissue\n' > "$issues"

while IFS= read -r folder || [[ -n "$folder" ]]; do
    folder=${folder%$'\r'}
    [[ -z "$folder" || "$folder" == \#* ]] && continue
    directory="$processed/$folder"
    info="$directory/contacts_unisex.info"
    fastq="$directory/R1.fq.gz"
    printf 'Processing %s\n' "$folder" >&2

    # One output row per folder; unavailable or invalid values are NA.
    fields=$'NA\tNA\tNA\tNA\tNA'
    if [[ -f "$info" && -r "$info" ]]; then
        if parsed=$(awk '
            NF {
                rows++
                if (NF != 5) bad=1
                for (i=2; i<=5; i++)
                    if ($i !~ /^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$/) bad=1
                value=sprintf("%s\t%s\t%s\t%s\t%s",$1,$2,$3,$4,$5)
            }
            END { if (rows != 1 || bad) exit 1; print value }
        ' "$info"); then
            fields=$parsed
        else
            printf '%s\tInvalid contacts_unisex.info: expected one nonblank row with five fields\n' "$folder" >> "$issues"
        fi
    else
        printf '%s\tMissing or unreadable contacts_unisex.info\n' "$folder" >> "$issues"
    fi

    read_count=NA
    if [[ -f "$fastq" && -r "$fastq" ]]; then
        # Standard FASTQ: four lines per read. No decompressed file is saved.
        # pipefail rejects gzip failures; awk rejects incomplete four-line groups.
        if count=$(gzip -cd -- "$fastq" | awk '
            END {
                if (NR % 4 != 0) exit 1
                printf "%.0f\n", NR / 4
            }
        '); then
            read_count=$count
        else
            printf '%s\tFASTQ decompression failed or line count was not divisible by four\n' "$folder" >> "$issues"
        fi
    else
        printf '%s\tMissing or unreadable R1.fq.gz\n' "$folder" >> "$issues"
    fi
    printf '%s\t%s\t%s\n' "$folder" "$fields" "$read_count" >> "$tmp"
done < "$folder_list"

mv -- "$issues" "${output%.tsv}.issues.tsv"
mv -- "$tmp" "$output"
printf 'Saved: %s\nIssues: %s\n' "$output" "${output%.tsv}.issues.tsv"
