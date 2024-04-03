#!/usr/bin/env python3

import json
import sys
import re


def get_novel_tx(class_code):
    """"""
    novel_class_codes = ["u", "y", "p"]

    if class_code[-2] in novel_class_codes:
        return True
    else:
        return False

def generate_xloc_2_gene_dict(loci_lines):
    """"""

    xloc_2_gene = {}            
    for i in range(len(loci_lines)):
        gene_col = loci_lines[i][2]
        all_ids = gene_col.split(",")
        gene_ids = []
        for id in all_ids:
            gene_ids.append(id.split("|")[0])
        uniq_ids = []
        for id in gene_ids:
            if id in uniq_ids:
                pass
            else:
                uniq_ids.append(id)

        xloc_2_gene[loci_lines[i][0]] = uniq_ids
    
    return xloc_2_gene

def generate_tx_ex_dict(gtf_list, sub_str, sub_str_1):
    """"""

    tmp_dict = {}
    tx_ids = []
    for i in range(len(gtf_list)):
        is_novel = False

        ##########################################
        # Searching for class code str in gtf data

        line   = gtf_list[i]
        fields = line.split("\t")

        tmp = re.compile(sub_str)
        hit = tmp.search(fields[8])

        if isinstance(hit, re.Match):
            class_code = str(hit.group(0)).strip()
        else:
            class_code = "class_code \"0\""
        

        is_novel = get_novel_tx(class_code)

        # Adding novel tx's
        if is_novel:
            tmp_dict[fields[8].split(sep=";")[1][-12:-1]] = [line]
            tx_ids.append(fields[8].split(";")[0])
        

        # Adding Exons to tx's
        elif not is_novel and isinstance(tmp_dict, dict) and fields[2] == 'exon' and fields[8].split(";")[0] in tx_ids:
            tmp  = re.compile(sub_str_1)
            hit  = tmp.search(line)
            xloc = str(hit.group(0))

            if xloc in tmp_dict.keys():
                tmp_value = tmp_dict[xloc]
                tmp_dict[xloc] = tmp_value + [line]
    
    return tmp_dict


def insert_novelty_tx(novel_tx_dict, loci_dict, overlap_tx_dict):
    """"""

    inserted_tx_ids = []
    for xloc, gtf_lines in novel_tx_dict.items():

            # If transcript xloc was aligned to some area of the genome
            if xloc in loci_dict.keys():

                # If there are transcripts that overlap only 1 reference gene and share an XLOC with our novel transcript
                if loci_dict[xloc][0] != "-" and len(loci_dict[xloc]) < 2 and loci_dict[xloc][0] in overlap_tx_dict.keys():
                    
                    for i in range(len(gtf_lines)):
                        fields       = gtf_lines[i].split("\t")
                        column_9     = fields[8].split(";")
                        if fields[2] == "transcript":
                            inserted_tx_ids.append(column_9[0])
                        tmp = column_9[0]
                        column_9[0] = column_9[1].strip()
                        column_9[1] = " " + tmp
                        column_9[0]  = "gene_id \"" + loci_dict[xloc][0] + "\""
                        fields[8]    = ";".join(column_9)
                        gtf_lines[i] = "\t".join(fields).strip()

                    old_values = overlap_tx_dict[loci_dict[xloc][0]]
                    overlap_tx_dict[loci_dict[xloc][0]] = old_values + gtf_lines
                    
      
    return overlap_tx_dict, inserted_tx_ids

def delete_inserted_lines(inserted_ids, gtf_lines):
    """"""

    new_gtf_lines = []
    for i in range(len(gtf_lines)):
        fields = gtf_lines[i].split("\t")
        column_9 = fields[8].split(";")
        if column_9[0] in inserted_ids:
            pass
        else:
            new_gtf_lines.append(gtf_lines[i])

    return new_gtf_lines

    

def main():

    #############################################################################
    ######################      Reading in work files     #######################

    try:
        combined_gtf = open(sys.argv[1], "r")
        gtf_list     = combined_gtf.readlines()

    except FileNotFoundError and IndexError:

        print("\nInput combined gtf file was not found or does not exist.\n")

    try:
        with open(sys.argv[2], "r") as updated_gene_id_file:
            updated_gene_id_gtf_dict = json.load(updated_gene_id_file)
            
    except FileNotFoundError and IndexError:

        print("\nInput gtf json file was not found or does not exist.\n")

    try:
        loci_file = open(sys.argv[3], "r")

        loci_lines = []
        for line in loci_file.readlines():
            loci_lines.append(line.split(sep="\t"))

    except FileNotFoundError and IndexError:

        print("\nInput gffcompare loci file was not found or does not exist.\n")

    try:
        gene_prefix = sys.argv[4]

    except IndexError:

        print("\nGene prefix identifier not passed into tool.\n")

    ##############################################################################
    ##############################################################################
        
    ################################################    
    #      Create dict with GENE prefix : GTF lines       #  
    
    sub_str   = ' class_code \".\"'
    sub_str_1 = gene_prefix + '_......'

    novel_tx_dict = generate_tx_ex_dict(gtf_list , sub_str, sub_str_1)
        
    # Dicts have been constructed for txs and exons #
    #################################################
        
    ###########################################
    #     Generate GENE prefix 2 GENE id Dict  &     #
    #   search updated gtf json for matches   #

    prefix_2_gene_dict = generate_xloc_2_gene_dict(loci_lines)

    #      GENE prefix 2 GENE ID dict generated      #
    ###########################################

    ####################################################
    # Insert Novel txs that share GENE prefix w/ overlap txs #

    inserted_novel_overlaps, deletion_ids = insert_novelty_tx(novel_tx_dict, prefix_2_gene_dict, updated_gene_id_gtf_dict)

    # Novel Txs that share GENE prefix w/ overlaps inserted  #
    ####################################################
    
    ###########################################################################
    # Delete all_double.combined.gtf lines that have been inserted into .json #
    
    new_gtf_list = delete_inserted_lines(deletion_ids, gtf_list)

    # Gtf List had been filtered for inserted novel transcripts to avoids duplication insert in later steps #
    #########################################################################################################

    ##########################################################################
    # Concatenating new gtf string and writing inserted .json + new gtf .gtf #

    gtf_str = ""
    for line in new_gtf_list:
        gtf_str = gtf_str + line

    with open("filtered_out_inserted_novels.gtf", "w") as new_gtf:
        new_gtf.write(gtf_str)
        new_gtf.close()

    with open("inserted_overlap_sharing_novels.json", "w") as new_json:
        json.dump(inserted_novel_overlaps, new_json, indent=4)
        


if __name__ == "__main__":

    main()











