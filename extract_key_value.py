#!/usr/bin/env python3 

import os
import json
import sys

path_to_tracking = sys.argv[1]

def extract_key_value(input_file):
    key_value_pairs = {}
    with open(input_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) >= 4:
                # Extract all parts between the second and the penultimate fields
                relevant_part = ' '.join(parts[2:-2])
                # Split by the first '|' symbol
                split_parts = relevant_part.split('|', 1)
                if len(split_parts) == 2:
                    value = split_parts[0].strip()
                    key = split_parts[1].strip()
                    key_value_pairs[key] = value

    return key_value_pairs

# Generate output JSON file name based on input file
output_json_name = os.path.basename(path_to_tracking).replace('.tracking', '_keyvalue.json')

extracted_data = extract_key_value(path_to_tracking)

with open("key_value_genes.json", 'w') as json_file:
    json.dump(extracted_data, json_file, indent=4)