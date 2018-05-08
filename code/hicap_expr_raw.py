#!/usr/bin/env python
# -*- coding: utf-8 -*

import os
import sys
import argparse


def int_dict_SP(int_files):

    """
        This function takes the interaction file of PD interaction and returns dictionary of gene name as
        the key and transcript_Id, distal chromosome position, as the list.

        Also it giives the dictionary of supporting pairs

        file = interaction file from HiCap runs with RefSeqName as the gene_ID, TranscriptName as transcript_ID

    """

    results_interaction_SP = {}

    for lines in int_files:
        line = lines.strip()
        fields = line.split("\t")

        if fields[11] == "-1":
            pass

        else:
            interaction_genes = fields[0]
            length = abs(int(fields[9])- int(fields[10]))
            interaction_status_SP = results_interaction_SP.get(interaction_genes,[])
            interaction_status_SP.append( [fields[12], fields[15],length])
            results_interaction_SP[interaction_genes] = interaction_status_SP


    return (results_interaction_SP)

def SP_CPM(results_interaction):

    """
        This is the extension of the above function that takes the dictinary of interction with their supporing
        pairs and retuns the counts per million of each interactor

    """
    counts_rep1 =0
    counts_rep2 =0
    final_list ={}

    for k in results_interaction.values():
        for value in k:
            counts_rep1 = int(counts_rep1) + int(value[0])
            counts_rep2 = int(counts_rep2) + int(value[1])

    for k in results_interaction.keys():
        Enh_rep1 = 0
        Enh_rep2 = 0
        tot_len_Enh = 0

        for values in results_interaction[k]:
            Enh_rep1 = int(values[0]) + Enh_rep1
            Enh_rep2 = int (values[1]) + Enh_rep2
            tot_len_Enh = int(values[2])+ tot_len_Enh

        tot_len_Enh = tot_len_Enh/1000
        Enh_rep1 = round(Enh_rep1/(counts_rep1/1000000),3)
        Enh_rep2 = round(Enh_rep2/(counts_rep2/1000000),3)
        #Enh_rep1 = round(Enh_rep1/((counts_rep1/1000000)*tot_len_Enh),3)
        #Enh_rep2 = round(Enh_rep2/((counts_rep2/1000000)*tot_len_Enh),3)
        final_list[k] = [Enh_rep1,Enh_rep2,tot_len_Enh]
    return (final_list, counts_rep1, counts_rep2)

def Main():
    parser = argparse.ArgumentParser(description="General command to manipulate interaction and expresstion data")
    parser.add_argument("PD", help = "PD dataset")
    parser.add_argument("-o", "--output", help ="output of interaction files", action='store', default=None)
    args = parser.parse_args()


    with open (args.PD , "r") as int_file_PD:
        next(int_file_PD)
        celltype_SP = int_dict_SP(int_file_PD)
        celltype_CPM_PD, replicate1_SP_count_PD, replicate2_SP_count_PD = SP_CPM(celltype_SP)

    print(replicate1_SP_count_PD,replicate2_SP_count_PD)


    if args.output:
        with open (args.output, "w") as out:
            for keys,values in celltype_CPM_PD.items():
                out.write (keys+"\t"+str(len(celltype_SP[keys]))+"\t"+"\t".join(str(i) for i in values))
                out.write ("\n")


if __name__ == "__main__":
    Main()