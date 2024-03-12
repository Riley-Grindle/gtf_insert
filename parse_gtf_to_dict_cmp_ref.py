#!/usr/bin/env python3 

import os
import json
import sys


path_to_gtf = sys.argv[1]

def parse_gtf_to_dict_by_cmp_ref(file_path):
    cmp_ref_data = {}
    with open(file_path, 'r') as file:
        current_cmp_ref = None
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            if fields[2] == "transcript" and 'cmp_ref' in fields[8]:
                cmp_ref_field = [f for f in fields[8].split(';') if 'cmp_ref' in f]
                if cmp_ref_field:
                    current_cmp_ref = cmp_ref_field[0].split('"')[1]
                    if current_cmp_ref not in cmp_ref_data:
                        cmp_ref_data[current_cmp_ref] = []

            if current_cmp_ref:
                cmp_ref_data[current_cmp_ref].append(line.strip())

    return cmp_ref_data

# Generate output JSON file name based on GTF file
output_json_name = os.path.basename(path_to_gtf).replace('.gtf', '_cmp_ref_data.json')

cmp_ref_dict = parse_gtf_to_dict_by_cmp_ref(path_to_gtf)

with open("gtf_to_dict_gffcmp.json", 'w') as json_file:
    json.dump(cmp_ref_dict, json_file, indent=4)