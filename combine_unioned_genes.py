#!/usr/bin/env python3

import json
import sys
import os
import re

def quicksort(arr):
    if len(arr) <= 1:
        return arr
    else:
        pivot = arr[0]
        left = [x for x in arr[1:] if int(x.split("\t")[3]) < int(pivot.split("\t")[3])]
        right = [x for x in arr[1:] if int(x.split("\t")[3]) >= int(pivot.split("\t")[3])]
        return quicksort(left) + [pivot] + quicksort(right)
 

def get_novel_tx(class_code):
    """"""
    novel_class_codes = ["r", "u", "i", "y", "p"]

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

def generate_tx_ex_dict(gtf_list, sub_str, sub_str_1):
    """"""

    tmp_dict = {}
    for i in range(len(gtf_list)):

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

        # Adding Exons to tx's
        elif not is_novel and isinstance(tmp_dict, dict) and fields[2] == 'exon':
            tmp  = re.compile(sub_str_1)
            hit  = tmp.search(line)
            xloc = str(hit.group(0))

            if xloc in tmp_dict.keys():
                tmp_value = tmp_dict[xloc]
                tmp_dict[xloc] = tmp_value + [line]
    
    return tmp_dict


def find_union_genes(loci_dict, gtf_dict):
    """"""

    # TO DO : 
    # --done-- identify values in loci_dict that are len() > 2
    # --done-- remove any keys that match those values above  
    # --done-- merge the values of removed keys under one key
    # --done-- the new key = the XLOC that is the loci_dict key
    # --done-- update the new gtf lines under the new XLOC key with gene_ids that match the XLOC, 
    # --done-- return unioned dict
    
    for key, value, in loci_dict.items():
        
        if len(value) > 1:
            if len(value) >= 6:
                print(key)
            tmp_value = []
            combined_ids = []
            for id in value:
                
                if id in gtf_dict.keys():
                    combined_ids.append(id)
                    tmp_value = tmp_value + gtf_dict[id]
                    del gtf_dict[id]

            for i in range(len(tmp_value)):
                line = tmp_value[i]
                fields = line.split("\t")
                column_9 = fields[8].split(";")
                counters = {}
                counter = 0
                for meta in column_9:
                    if "gene_id" in meta:
                        meta = "gene_id \"" + key + "\""
                        counters[counter] = meta
                    if " gene_name" in meta:
                        meta = " gene_name \"" + key + "\""
                        counters[counter] = meta
                    counter += 1
                for index in list(counters.keys()):
                    column_9[index] = counters[index]
                fields[8] = ";".join(column_9)
                line = "\t".join(fields).strip()
                tmp_value[i] = line
            
            combined_str = ''
            for id in combined_ids:
                combined_str = combined_str + id + ","

            gtf_dict[key] = [tmp_value, combined_str[:-1]]
    
    
    
    return gtf_dict


def insert_novelty_tx(novel_tx_dict, combined_tx_dict):
    """"""

    # TO DO : 
    # ------- identify XLOC keys in novel_tx_dict that are conserved in combined_tx_dict
    # ------- insert gtf lines from novel dict into value of combined dict
    # ------- add transcript id to list of ids that need to be deleted from gffcmp gtf
    # ------- return inserted dict and inserted tx ids

    inserted_ids = []
    for xloc, gtf_lines in novel_tx_dict.items():

        if xloc in combined_tx_dict.keys():

            old_values = combined_tx_dict[xloc]
            new_values = []
            for line in gtf_lines:
                fields = line.split("\t")
                column_9 = fields[8].split(";")
                tmp = column_9[0]
                column_9[0] = column_9[1].strip()
                column_9[1] = " " + tmp
                fields[8] = ";".join(column_9)
                line = "\t".join(fields)
                new_values.append(line)

            combined_tx_dict[xloc] = old_values + new_values

            for i in range(len(gtf_lines)):
                fields       = gtf_lines[i].split("\t")
                column_9     = fields[8].split(";")
                if fields[2] == "transcript":
                    inserted_ids.append(column_9[0])

    
    return combined_tx_dict, inserted_ids

def sort_xlocs(fully_combined_inserted_dict, prefix):

    # TO DO : 
    # ------- Identify keys that are XLOCs and sort values at the transcript level
    # ------- re-insert exons in ascending order using tx id 
    # ------- delete gene level annotation and insert new gene annotation covering all txs
    # ------- add parent gene data showing which parent genes contributed to union gene
    # ------- return final_sorted_dict

    counter = 0
    for key, value in fully_combined_inserted_dict.items():

        if prefix in key:
            counter += 1

            sorted_line_list = []
            other_list = []
            for i in range(len(value[0])):

                line = value[0][i]
                fields = line.split("\t")
                
                if fields[2] == "gene":
                    pass
                elif fields[2] in ['transcript', "RNA", "mRNA"]:
                    
                    if len(sorted_line_list) > 0:
                        sorted_line_list.append(line.strip())

                    else:
                        new_gene_level = fields[0] + "\tStringTie" + "\tgene" + "\t0" + "\t0" + "\t" + fields[5] + "\t" + fields[6] + "\t" + fields[7] + "\tgene_id " + "\"" + key + "\"" + "; gene_name " + "\"" + key + "\"" + "; combined_IDs " + "\"" + value[1] + "\""
                        sorted_line_list.append(new_gene_level)
                        sorted_line_list.append(line.strip())
                else:
                    other_list.append(line)
            
            sorted_line_list = quicksort(sorted_line_list)
            gene_fields = sorted_line_list[0].split("\t")
            gene_fields[3] = sorted_line_list[1].split("\t")[3]
            max_boundary = 0
            for i in range(len(sorted_line_list)):
                if int(sorted_line_list[i].split("\t")[4]) > max_boundary:
                    max_boundary = int(sorted_line_list[i].split("\t")[4])
            gene_fields[4] = str(max_boundary)
            sorted_line_list[0] = "\t".join(gene_fields)    
            
            for line in other_list[::-1]:
                fields = line.split("\t")
                column_9 = fields[8].split(";")
                tx_index = [index for index, item in enumerate(column_9) if "transcript_id" in item]
                tx_id = column_9[tx_index[0]]

                index = 0
                for tx_line in sorted_line_list:

                    if tx_id in tx_line:
                        if index == 0:
                            sorted_line_list = [sorted_line_list[0]] + [line.strip()] + sorted_line_list[1:]
                            #print(line)
                            #print(sorted_line_list)
                            break
                        else:
                            sorted_line_list = sorted_line_list[:index+1] + [line.strip()] + sorted_line_list[index+1:]
                            #print(line)
                            #print(sorted_line_list)
                            break

                    index += 1

            fully_combined_inserted_dict[key] = sorted_line_list


    return fully_combined_inserted_dict


def main():
    
    #############################################################################
    ######################      Reading in work files     #######################

    try:
        combined_gtf = open(sys.argv[1], "r")
        gtf_list     = combined_gtf.readlines()

    except FileNotFoundError and IndexError:

        print("\nInput combined gtf file was not found or does not exist.\n")

    try:
        with open(sys.argv[2], "r") as reference_gtf_inserted_file:
            reference_inserted_dict = json.load(reference_gtf_inserted_file)
            
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
    sub_str_1 = gene_prefix + 'G_......'

    novel_tx_dict = generate_tx_ex_dict(gtf_list , sub_str, sub_str_1)
        
    # Dicts have been constructed for txs and exons #
    #################################################

    ###########################################
    #     Generate GENE prefix 2 GENE id Dict  &     #
    #   search updated gtf json for matches   #

    
    gene_2_gene_dict = generate_xloc_2_gene_dict(loci_lines)

    #      GENE prefix 2 GENE ID dict generated      #
    ###########################################

    #########################################################################
    # Update gtf genes to GENE prefix if more than one gene is spanned by an GENE prefix #

    unionized_dict  = find_union_genes(gene_2_gene_dict, reference_inserted_dict)

    #                                                                        #
    ##########################################################################

    ################################################################
    # Insert novel transcripts that share GENE prefix with unioned genes #

    inserted_novel_tx_dict, deletion_ids = insert_novelty_tx(novel_tx_dict, unionized_dict)

    #                                                              #
    ################################################################

    ###########################################################################################################
    # Sort unionized genes, remove uneeded gene levels, add total bounding gene level with updated start-stop #

    sorted_inserted_novels_overlaps_refs = sort_xlocs(inserted_novel_tx_dict, gene_prefix)

    ###########################################################################
    # Delete all_double.combined.gtf lines that have been inserted into .json #
    
    new_gtf_list = delete_inserted_lines(deletion_ids, gtf_list)

    # Gtf List had been filtered for inserted novel transcripts to avoids duplication insert in later steps #
    #########################################################################################################

    
    gtf_str = ""
    for line in new_gtf_list:
        gtf_str = gtf_str + line

    with open("filtered_out_inserted_novels_2.gtf", "w") as new_gtf:
        new_gtf.write(gtf_str)
        new_gtf.close()

    with open("joined_genes_w_spanning_txs.json", "w") as new_json:
        json.dump(sorted_inserted_novels_overlaps_refs, new_json, indent=4)
    

if __name__ == "__main__":

    main()
