#!/usr/bin/env python3

import os
import json
import sys

path_to_gtf = sys.argv[1]


# Assuming overlapClassCodes is passed as '["k", "m", "n", "j", "e"]'
overlap_class_codes = ["k", "m", "n", "j", "e", "o"]


def find_overlap_transcript(gtf, overlap_class_codes):
    # Generate output file name
    output_gtf_name = os.path.basename(gtf).replace('.gtf', '_ready_insert_transcript.gtf')

    # Open the input GTF file
    with open(gtf, 'r') as file:
        combined_gtf_lines = file.readlines()

    overlap_transcript_lines = []
    for line in combined_gtf_lines:
        fields = line.strip().split('\t')
        if len(fields) < 9:
            continue

        if fields[2] == "transcript" and 'class_code "' in fields[8]:
            class_code = fields[8].split('class_code "')[1][0]
            if class_code in overlap_class_codes:
                overlap_transcript_lines.append(line)
                # Add all subsequent exon lines
                index = combined_gtf_lines.index(line) + 1
                while index < len(combined_gtf_lines) and 'exon' in combined_gtf_lines[index].split('\t')[2]:
                    overlap_transcript_lines.append(combined_gtf_lines[index])
                    index += 1

    # Write the overlap transcript lines to the output GTF file
    with open("find_overlaps.gtf", 'w') as output_file:
        output_file.writelines(overlap_transcript_lines)

    return output_gtf_name

# Call the function and pass the list of overlap class codes
output_gtf_path = find_overlap_transcript(path_to_gtf, overlap_class_codes)
