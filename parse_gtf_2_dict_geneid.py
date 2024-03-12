#!/usr/bin/env python3 

import os
import json
import sys

path_to_gtf = sys.argv[1]


def parse_gtf_to_dict_by_geneID(file_path):
    gene_data = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue  # Skip comment lines
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            gene_id_field = [f for f in fields[8].split(';') if 'gene_id' in f]
            if gene_id_field:
                try:
                    gene_id = gene_id_field[0].split('"')[1]
                except IndexError:
                    gene_id = gene_id_field[0].split('=')[1]

                if gene_id not in gene_data:
                    gene_data[gene_id] = []

                gene_data[gene_id].append(line.strip())

    return gene_data

# Generate output JSON file name based on GTF file
output_json_name = os.path.basename(path_to_gtf).replace('.gtf', '_gene_id_data.json')

gene_id_dict = parse_gtf_to_dict_by_geneID(path_to_gtf)

with open("parse_reference_gtf.json", 'w') as json_file:
    json.dump(gene_id_dict, json_file, indent=4)