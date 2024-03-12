#!/usr/bin/env python3 

import os
import json
import sys

path_to_key_values = sys.argv[1]
path_to_overlap_json = sys.argv[2]

def update_keys(source_file, target_file):
    try:
        with open(source_file, 'r') as file:
            source_data = json.load(file)

        with open(target_file, 'r') as file:
            target_data = json.load(file)

        updated_data = {}
        for key, value in target_data.items():
            if key in source_data:
                # If the key exists in source_data, update the key
                new_key = source_data[key]
                new_value = []
                for item in value:
                    field_list = item.split(sep = "\t")
                    column_nine = field_list[8].split(";")
                    column_nine[1] = "gene_id " + "\"" + source_data[key] + "\""
                    tmp = column_nine[0]
                    column_nine[0] = column_nine[1].strip()
                    column_nine[1] = " " + tmp
                    column_nine = ";".join(column_nine)
                    field_list[8] = column_nine
                    fields = "\t".join(field_list)
                    new_value.append(fields)
                # Use setdefault to append values for the same new_key
                updated_data.setdefault(new_key, []).extend(new_value)

        # Generate output JSON file name based on target file
        output_json_name = os.path.basename(target_file).replace('.json', '_updated.json')
        with open("updated_gffcmp.json", 'w') as file:
            json.dump(updated_data, file, indent=4)

        return "Update successful and saved to new file"
    except Exception as e:
        return f"Error: {e}"

result = update_keys(path_to_key_values, path_to_overlap_json)
print(result)
