#!/usr/bin/python
#! -*- coding:utf-8 -*-

import argparse
def converInteractionTable(filename,savefilename):
    with open(filename,'r') as file_handle:
        lines = file_handle.readlines()

    save_file = open(savefilename,'w')
    for line in lines[1:]:
        line = line[:len(line)-1]
        elements = line.split(';')
        bacteria_id = elements[0]
        bacterial_name = elements[1]
        specie = elements[2]
        genus = elements[3]
        phage_id  = elements[4]
        phage_name = elements[5]
        spacer_id = elements[6]
        missmatch = elements[7]
        phage_position = elements[8]
        flag = elements[9]
        ## 填写sql 的Insert语句
        sql = "INSERT INTO `interaction` \
(`bacteria_id`, `bacteria_name`, `specie`, \
`genus`, `phage_id`, `phage_name`, `spacer_id`, \
`missmatch`, `phage_position`, `flag`)\
VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" %(\
            bacteria_id,bacterial_name,specie,genus,phage_id,phage_name,spacer_id,missmatch,phage_position,flag )
        save_file.write(sql+'\n')
    save_file.close()

def converHitProteinTable(filename,savefilename):
    with open(filename,'r') as file_handle:
        lines = file_handle.readlines()

    save_file = open(savefilename,'w')
    for line in lines[1:]:
        line = line[:len(line)-1]
        elements = line.split(';')
        print(line,elements)
        bac_id = elements[0]
        phage_id = elements[1]
        spacer_id = elements[2]
        protein_id = elements[3]
        protein_position = elements[4]
        protein_product = elements[5]
        protein_sequence = elements[6]

        insertsql = "INSERT INTO `hitprotein`\
 (`bac_id`, `phage_id`, `spacer_id`, `protein_id`, \
 `protein_position`, `protein_product`, `protein_sequence`)\
VALUES ('%s','%s','%s','%s','%s','%s','%s')" %(bac_id,phage_id,spacer_id,protein_id,protein_position,protein_product,protein_sequence)
        save_file.write(insertsql+'\n')
    save_file.close()
def converSequencesTable(filename,savefilename):
    with open(filename,'r') as file_handle:
        lines = file_handle.readlines()

    save_file = open(savefilename,'w')
    for line in lines[1:]:
        line = line[:len(line)-1]
        elements = line.split(';')
        bac_id = elements[0]
        phage_id = elements[1]
        spacer_id = elements[2]
        spacer_sequence = elements[3]
        hit_sequence = elements[4]
        long_hit_sequence = elements[5]
        coverage = elements[6]

        insertsql = "INSERT INTO `sequence` (`bac_id`, `phage_id`, `spacer_id`, `spacer_sequence`, `hit_sequence`, `long_hit_sequence`, `coverage`)\
VALUES ('%s','%s','%s','%s','%s','%s','%s') " %(bac_id,phage_id,spacer_id,spacer_sequence,hit_sequence,long_hit_sequence,coverage)
        save_file.write(insertsql+'\n')
    save_file.close()
def init():
    parser = argparse.ArgumentParser(description='conver the csv to sql table')

    parser.add_argument('-i', '--input', help="the csv file")
    parser.add_argument('-s', '--savefile', help="the sql file to save")
    parser.add_argument('-t', '--type', type = str,choices=['sequences', 'hitProtein', 'interaction'],)
    args = parser.parse_args()
    if args.input:
        inputfile = args.input
    else:
        print("please input the right args")
        exit(0)
    if args.savefile:
        savefile = args.savefile
    else:
        print("please input the right args")
        exit(0)
    if args.type:
        type = args.type
    else:
        print("please input the right args")
        exit(0)
    return (inputfile,savefile,type)
if __name__=="__main__":

    inputfile,savefile,type =  init()
    if type=="sequences":
        converSequencesTable(inputfile,savefile)
    if type=="hitProtein":
        converHitProteinTable(inputfile,savefile)
    if type =="interaction":
        converInteractionTable(inputfile,savefile)
