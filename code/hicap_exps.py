#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse


"""
    This is a general functions developed for manipulating interaction data and expression data

     python hicap_exps.py /Volumes/Work_drive/prj/Expression_Vs_Interaction/data/THP1/gene_expression_THP1_nlps_filtered.txt /Volumes/Work_drive/prj/Expression_Vs_Interaction/data/THP1/THP1run.hg19.nLPS.Proximities.Probe_Distal.filtered.txt /Volumes/Work_drive/prj/Expression_Vs_Interaction/data/THP1/THP1run.hg19.nLPS.Proximities.Probe_Probe.SMCscaled.filtered.txt

"""


def exp_dict(exp_data):

    """
        The function take expression file given input and returns dictionary of the expression data as the values for
        every genes. At the mone these are genes that are expressed in either of the replicates
        Since the expression values are given here in the last two columns, I have taken these last two columns

    """

    dict_exp = {}
    next(exp_data)

    for lines in exp_data:
        line = lines.strip().split("\t")
        dict_fields = line[0]
        fields = dict_exp.get(dict_fields,"")
        fields = line[-2],line[-1]
        dict_exp[dict_fields] = fields
    return(dict_exp)

def int_dict_pd(int_files):

    """
        This function takes the interaction file of PD interaction and returns dictionary of gene name as
        the key and transcript_Id, distal chromosome position, as the list.

        Also it giives the dictionary of supporting pairs

        file = interaction file from HiCap runs with RefSeqName as the gene_ID, TranscriptName as transcript_ID

    """

    results_interaction = {}
    results_interaction_SP = {}

    for lines in int_files:
        line = lines.strip()
        fields = line.split("\t")

        if fields[11] == "-1":
            pass

        else:
            interaction_genes = fields[0]
            interaction_status_genes = results_interaction.get(interaction_genes,[])

            interaction_status_genes.append([fields[1], fields[8], fields[9], fields[10], fields[11]])
            results_interaction[interaction_genes] = interaction_status_genes

            length = abs(int(fields[9])- int(fields[10]))
            interaction_status_SP = results_interaction_SP.get(interaction_genes,[])
            interaction_status_SP.append( [fields[12], fields[15],length])
            results_interaction_SP[interaction_genes] = interaction_status_SP


    return (results_interaction, results_interaction_SP)

def int_dict_pp(int_files_PP):

    """
        This function takes the interaction file of PP interaction and returns dictionary of gene name as
        the key and transcript_Id, distal chromosome position, as the list.

        file = interaction file from HiCap runs with RefSeqName as the gene_ID, TranscriptName as transcript_ID

    """
    results_interact_PP = {}
    results_interact_PP_SP = {}

    for lines in int_files_PP:
        line = lines.strip()
        fields = line.split("\t")
        interaction_genes = fields[0]

        interaction_status_genes = results_interact_PP.get(interaction_genes,[])
        interaction_status_genes.append([fields[1], fields[11], fields[12], fields[13], fields[16]])
        results_interact_PP[interaction_genes] = interaction_status_genes

        length = abs(int(fields[12])- int(fields[13]))
        interaction_status_SP = results_interact_PP_SP.get(interaction_genes,[])
        interaction_status_SP.append([fields[17], fields[20],length])
        results_interact_PP_SP[interaction_genes] = interaction_status_SP

    return (results_interact_PP, results_interact_PP_SP)

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
    parser.add_argument("exprs", help = "Expression file of given celltype; Currently takes only expressed dataset")
    parser.add_argument("PD", help = "PD dataset")
    parser.add_argument("PP", help = "PP dataset")

#    parser.add_argument("-o", "--output", help ="output of interaction files", action='store', default=None)

    args = parser.parse_args()

    with open (args.exprs, "r") as exp_fil:
        celltype_expr_dict = exp_dict(exp_fil)

    with open (args.PD , "r") as int_file_PD:
        next(int_file_PD)
        celltype_pd, celltype_SP = int_dict_pd(int_file_PD)
        celltype_CPM_PD, replicate1_SP_count_PD, replicate2_SP_count_PD = SP_CPM(celltype_SP)

    print(replicate1_SP_count_PD,replicate2_SP_count_PD)


    with open (args.PP , "r") as int_file_PP:
        next(int_file_PP)
        celltype_pp, celltype_pp_SP = int_dict_pp(int_file_PP)
        celltype_CPM_PP, replicate1_SP_count_PP, replicate2_SP_count_PP = SP_CPM(celltype_pp_SP)

    print(replicate1_SP_count_PP,replicate2_SP_count_PP)


if __name__ == "__main__":
    Main()