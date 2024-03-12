#!/usr/bin/env python3 

import os
import json
import sys

path_to_new_json = sys.argv[1]
path_to_ref_json = sys.argv[2]


def separate_transcripts_with_same_cmp_ref_geneID(gtf_list):
    processed_dict = {}
    key_counter = 1
    key = 1
    for line in gtf_list:
        fields = line.split('\t')
        if fields[2] in ['transcript', 'exon', 'mRNA', 'RNA']:
            if fields[2] == 'transcript' or fields[2] == 'mRNA' or fields[2] == 'RNA':
                start_position = fields[3]
            else:
                start_position
            #start_position = fields[3] if fields[2] == 'transcript' or fields[2] == 'mRNA' else start_position
            processed_dict.setdefault(start_position, []).append(line)
        elif fields[2] == 'gene':
            processed_dict['gene'] = [line]
        else:
            unique_key = f'unique_key_{key_counter}'
            processed_dict[unique_key] = [line]
            key_counter += 1

    return processed_dict

def find_last_transcript_key(reference_value):
    last_transcript_key = None
    for key in reference_value:
        if key.startswith("unique_key_"):
            continue
        last_transcript_key = key
    return last_transcript_key

# Load JSON files
with open(path_to_new_json, 'r') as file:
    cmp_ref_data = json.load(file)

with open(path_to_ref_json, 'r') as file:
    reference_data = json.load(file)

# Process each cmp_ref_data
count = 0
for cmp_ref, lines in cmp_ref_data.items():
    
    if cmp_ref in reference_data:
        
        cmp_ref_value = separate_transcripts_with_same_cmp_ref_geneID(lines)
        reference_value = separate_transcripts_with_same_cmp_ref_geneID(reference_data[cmp_ref])
        
        for cmp_ref_start_position, cmp_ref_lines in sorted(cmp_ref_value.items(), key=lambda item: item[0], reverse=True):
            inserted = False
            for reference_start_position, reference_lines in reference_value.items():
                
                # Updating gene start/stop if transcript extends it
                if reference_start_position == 'gene':

                    fields = reference_lines[0].split('\t')
                    if int(fields[3]) >= int(cmp_ref_start_position):
                        fields[3] = cmp_ref_start_position
                        start_value = cmp_ref_start_position
                        

                    new_fields = cmp_ref_lines[-1].split('\t')
                    if int(fields[4]) <= int(new_fields[4]):
                        fields[4] = new_fields[4]
                        stop_value = new_fields[4]
                        

                    reference_lines[0] = "\t".join(fields)
                    reference_value[reference_start_position] = reference_lines

                

                if reference_start_position.isdigit() and int(cmp_ref_start_position) <= int(reference_start_position):
                    reference_value[reference_start_position] = cmp_ref_lines + reference_lines
                    inserted = True
                    break
                

            if not inserted:
                last_transcript_key = find_last_transcript_key(reference_value)
                if last_transcript_key is not None:
                    reference_value[last_transcript_key].extend(cmp_ref_lines)
                else:
                    unique_key = f'unique_end_key_{cmp_ref_start_position}'
                    reference_value[unique_key] = cmp_ref_lines


        reference_data[cmp_ref] = [line for sublist in reference_value.values() for line in sublist]
    



# Generate output JSON file name
output_json_name = 'finished_insert_overlap_transcript.json'
with open("overlap_inserted.json", 'w') as json_file:
    json.dump(reference_data, json_file, indent=4)

