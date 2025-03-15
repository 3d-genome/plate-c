# obtain list of protein-coding genes from 10x reference
cd ~/tools/cellranger/references/refdata-gex-GRCm39-2024-A/genes
zcat genes.gtf.gz | awk -v FS='\t' '($3=="gene")' | grep "gene_type \"protein_coding\";" | grep -o 'gene_name "[^"]*"' | cut -d '"' -f 2 | sort | uniq > refdata-gex-GRCm39-2024-A.gene_names.protein_coding.txt
