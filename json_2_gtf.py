#!/usr/bin/env python3 

import os
import json
import sys

path_to_updated_overlaps_json = sys.argv[1]

def json_to_gtf(json_path, output_gtf_path):
    with open(json_path, 'r') as json_file:
        gene_data = json.load(json_file)

    with open(output_gtf_path, 'w') as gtf_file:
        for gene_id in gene_data:
            for line in gene_data[gene_id]:
                gtf_file.write(line + '\n')

# Generate output GTF file name based on JSON file
output_gtf_name ='overlap_inserted.gtf'

json_to_gtf(path_to_updated_overlaps_json, output_gtf_name)