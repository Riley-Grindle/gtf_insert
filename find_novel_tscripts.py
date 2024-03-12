#!/usr/bin/env python3 

import os
import json
import sys

path_to_gffcmp_gtf = sys.argv[1]
path_to_newest_gtf = sys.argv[2]
path_to_tracking = sys.argv[3]

# Convert the JSON string to a Python list
novel_class_codes = ["r", "u", "i", "y", "p"]

def find_novel_transcripts(gtf, reference_gtf_path, novel_class_codes, tracking):


    # Generate output file name based on combined GTF file
    output_gtf_name ='final_annotation.gtf'

    with open(tracking, "r") as file:
        tracking_dict = json.load(file)
    tx_keys = list(tracking_dict.keys())

    # Read the combined GTF file
    with open(gtf, 'r') as file:
        combined_gtf_lines = file.readlines()

    novel_transcript_lines = []
    target_class_code = False
    for line in combined_gtf_lines:
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue

        # Check for transcript lines with novel class codes
        if fields[2] == "transcript" and 'class_code "' in fields[8]:
            class_code = fields[8].split('class_code "')[1][0]
            target_class_code = class_code in novel_class_codes
            if target_class_code:
                column_nine = fields[8].split(";")
                if 'cmp_ref "' in fields[8]:
                    tmp = column_nine[3]
                    if tmp[10:-1] in tx_keys:
                         column_nine[1] = "gene_id " + "\"" + tracking_dict[column_nine[3][10:-1]] + "\""
                         tmp = column_nine[0]
                         column_nine[0] = column_nine[1].strip()
                         column_nine[1] = " " + tmp
                    column_nine = ";".join(column_nine)
                    fields[8] = column_nine
                else:
                    tmp = column_nine[0]
                    column_nine[0] = column_nine[1].strip()
                    column_nine[1] = " " + tmp
                    column_nine = ";".join(column_nine)
                    fields[8] = column_nine
                new_line = "\t".join(fields) + "\n"
                novel_transcript_lines.append(new_line)

        elif target_class_code:
            column_nine = fields[8].split(";")
            tmp = column_nine[0]
            column_nine[0] = column_nine[1].strip()
            column_nine[1] = " " + tmp
            column_nine = ";".join(column_nine)
            fields[8] = column_nine
            new_line = "\t".join(fields) + "\n"
            novel_transcript_lines.append(new_line)

    # Read the reference GTF file and append the novel transcript lines
    with open(reference_gtf_path, 'r') as ref_file:
        reference_gtf_lines = ref_file.readlines()

    # Write the combined reference and novel transcript lines to the output GTF file
    with open("final_annotation.gtf", 'w') as output_file:
        output_file.writelines(reference_gtf_lines + novel_transcript_lines)

    return output_gtf_name

output_gtf_path = find_novel_transcripts(path_to_gffcmp_gtf, path_to_newest_gtf,novel_class_codes, path_to_tracking)
