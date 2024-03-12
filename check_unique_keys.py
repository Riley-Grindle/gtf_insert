#!/usr/bin/env python3 

import os
import json
import sys

path_to_transcript_gene_map = sys.argv[1]

def check_unique_keys(json_file):
    try:
        with open(json_file, 'r') as file:
            data = json.load(file)

        keys = list(data.keys())
        if len(keys) == len(set(keys)):
            return True
        else:
            return False
    except Exception as e:
        return str(e)

result = check_unique_keys(path_to_transcript_gene_map)
result_text = f"Are all keys unique? {result}"

# Save the result to a text file
output_text_name = 'results.txt'
with open(output_text_name, 'w') as output_file:
    output_file.write(result_text)