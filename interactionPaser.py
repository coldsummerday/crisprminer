#! /usr/bin/env python
# ! -*- code:utf-8 -*-
from Bio import Seq
import os
from Bio import SeqIO
from Bio import Entrez
from time import sleep
import time

Entrez.email = "coldsummerday1211@gmail.com"

phageInfo_cache = {}


def get_file_path_hash(refseqid):
    number = 0
    for a in refseqid:
        number += ord(a)
    saveNumber = number % 150
    return str(saveNumber)


def get_spc_sequence(spcfile, spacerid):
    spacerid = ">" + spacerid
    with open(spcfile, 'r') as  file_handle:
        lines = file_handle.readline()
    for index, line in enumerate(lines):
        line = line[0:len(line) - 1]
        if line == spacerid:
            seq = lines[index + 1]
            if seq[-1] == '\n':
                seq = seq[0:len(seq) - 1]
        return seq


def get_genus_specis(speciefile):
    with open(speciefile, 'r') as file_handle:
        lines = file_handle.readlines()
    bacteriaToGenus = {}
    for line in lines:
        if line[-1] == '\n' or line[-1] == '\t':
            line = line[:len(line) - 1]
        elements = line.split(";")
        id = elements[3]
        genus = elements[1]
        specie = elements[2]
        bacteriaToGenus[id] = {'genus': genus, 'specie': specie}
    return bacteriaToGenus


def get_id_to_name(namefile):
    with open(namefile) as file_handle:
        lines = file_handle.readlines()
    bacidToName = {}
    for line in lines:
        if line[-1] == '\n' or line[-1] == '\t':
            line = line[:len(line) - 1]
        element = line.split("---- ")
        bacid = element[0]
        name = element[1]
        bacidToName[bacid] = name


class BacterialAsmParser(object):
    def __init__(self, bacterialId, asmSavePath):
        self.bacterialPath = '/zrom/jobs/bacteria/'
        self.bacterialId = bacterialId
        self.asmPath = asmSavePath
        self.interaction_list = []
        self.parserFile()
        self.genus = None
        self.specie = None
        self.get_bac_name()
        self.get_specie_genus()

    def get_specie_genus(self):
        global idToGenus_specisDic
        if self.bacterialId in idToGenus_specisDic.keys():
            self.genus = idToGenus_specisDic[id]['genus']
            self.specie = idToGenus_specisDic[id]['specie']

    def get_bac_name(self):
        global idToNameDic
        if self.bacterialId in idToNameDic.keys():
            self.bacterialName = idToNameDic[self.bacterialId]

    def parserFile(self):
        with open(self.asmPath + self.bacterialId + '.asn') as file_handle:
            lines = file_handle.readlines()
        if len(lines) == 0:
            self.empty_flag = True
            return
        for line in lines:
            elements = line.split('\t')
            spacerid = elements[0]
            phageid = elements[1].split('.')[0]
            identy = float(elements[2])
            alength = float(elements[3])
            missmatch = float(elements[4])
            start = int(elements[8])
            end = int(elements[9])
            length = alength / identy
            if length < 22:
                self.interaction_list.append \
                    ({'spacerid': spacerid, 'phageid': phageid, \
                      'missmatch': missmatch, 'start': start, 'end': end})
                continue
            if missmatch < ((length - 22) ** 0.5):
                self.interaction_list.append \
                    ({'spacerid': spacerid, 'phageid': phageid, \
                      'missmatch': missmatch, 'start': start, 'end': end})


class phageInfo(object):
    def __init__(self, phageID, phagefilePath):
        self.phageID = phageID
        self.phagefilePath = phagefilePath
        self.proteinList = []
        self.phageName = None
        self.phageSeq = ''
        self.get_phage_seq()
        self.parserGenbankFile()

    def get_phage_seq(self):
        with open(self.phagefilePath + self.phageID + '.fna', 'r') as file_handle:
            lines = file_handle.readlines()
        alllines = ''
        for line in lines:
            if line[0] == '>':
                continue
            line = line[0:len(line) - 1]
            alllines += line
        self.phageSeq = alllines

    def get_position_seq(self, start, end):
        return self.phageSeq[start:end]

    def parserGenbankFile(self):
        if not os.path.exists(self.phagefilePath + self.phageID + '.gb'):
            try:
                sleep(0.15)
                handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=self.phageID)
                with open(self.phagefilePath + self.phageID + '.gb', 'w') as file_handle:
                    file_handle.write(handle.read())
            except:
                return
        records = list(SeqIO.parse(self.phagefilePath + self.phageID + '.gb', "genbank"))
        for record in records:
            self.phageName = record.annotations['organism']
            for feature in record.features:
                if feature.type == "CDS":
                    location = str(feature.location)
                    product = feature.qualifiers['product'][0] \
                        if feature.qualifiers.has_key('product') else 'unkown'
                    proteinId = feature.qualifiers['protein_id'][0] \
                        if feature.qualifiers.has_key(
                        'protein_id') else 'unkown'
                    translation = feature.qualifiers['translation'][0] \
                        if feature.qualifiers.has_key(
                        'translation') else 'unkown'
                    protein = phageProtein(proteinId, location, product, translation)
                    protein.parserLocation()
                    self.proteinList.append(protein)

    def get_interaction_protein(self, start, end, order):
        for protein in self.proteinList:
            if order == '+' and protein.order == '+':
                if start > protein.start and end < protein.end:
                    return protein
            elif order == '-' and protein.order == '-':
                if start > protein.start and end < protein.end:
                    return protein
        return None


class phageProtein(object):
    def __init__(self, proteinID, location, product, seq):
        self.proteinID = proteinID
        self.location = location
        self.product = product
        self.seq = seq

    def parserLocation(self):
        start = int(self.location.split(':')[0].split('[')[1])
        end = int(self.location.split(':')[1].split(']')[0])
        self.order = self.location.split('(')[1].split(')')[0]
        if start > end:
            temp = start
            start = end
            end = temp
        self.start = start
        self.end = end


def get_phageInfo_from_cache(phageid):
    global phageInfo_cache
    global phage_file_path
    if phageid in phageInfo_cache.keys():
        return phageInfo_cache[phageid]
    else:
        phage_handle = phageInfo(phageid, phage_file_path + phageid + '/')
        phageInfo_cache[phageid] = phage_handle
        return phage_handle


def paserAsn(asnFileDir):
    global hitproteinfile_handle
    global seqfile_handle
    global interactionsfile_handle
    filenamelist = os.listdir(asnFileDir)
    filenumber = len(filenamelist)
    for index, filename in enumerate(filenamelist):
        if os.path.getsize(asnFileDir + filename) == 0:
            continue
        bacid = filename.split('.')[0]
        parser_handle = BacterialAsmParser(bacid, asnFileDir)
        for interaction in parser_handle.interaction_list:
            start = interaction['start']
            end = interaction['end']
            spacer_id = interaction['spacerid']
            missmatch = interaction['missmatch']
            phageid = interaction['phageid']
            phageInfo_handle = get_phageInfo_from_cache(phageid)
            if start < end:
                order = '+'
                hitsequens = phageInfo_handle.get_position_seq(int(start) - 1,int(end))

                longSequens = phageInfo_handle.get_position_seq(int(start) - 6,int(end) + 5)
            else:
                temp = start
                start = end
                end = temp
                order = '-'
                sequens = str(Seq(phageInfo_handle.get_position_seq(int(start) - 1,int(end))).reverse_complement())
                longSequens = str(Seq(phageInfo_handle.get_position_seq(int(start) - 6,int(end) + 5)).reverse_complement())
            spacer_seq = get_spc_sequence(parser_handle.bacterialPath+get_file_path_hash(bacid) + \
                                          '/' + bacid + '/' + bacid + '.spc',spacer_id )
            seqfile_handle.write(bacid+';'+ str(phageInfo_handle.phageID) + ';' + \
                                spacer_id + ';' + str(spacer_seq) + ';' + str(hitsequens) + \
                                 ";" + longSequens + ';' + str(float(len(hitsequens))/float(len(spacer_seq))) + '\n')
            seqfile_handle.flush()
            hitprotein = phageInfo_handle.get_interaction_protein(start, end, order)
            if not hitprotein:
                flag = 0
            else:
                # 代表匹配那段基因有编码的蛋白质
                flag = 1
                hitproteinfile_handle.write(bacid + ';' + phageid + ';' + \
                                            str(hitprotein.proteinID) + ';' + \
                                            str(hitprotein.location) + ';' + \
                                            str(hitprotein.product) + ';' + \
                                            str(hitprotein.seq) + '\n')
                hitproteinfile_handle.flush()
            interactionsfile_handle.write(bacid + \
                                          ';' + str(parser_handle.bacterialName) + \
                                          str(parser_handle.specie) + ';' + \
                                          str(parser_handle.genus) + ';' + \
                                          str(phageid) + ';' + \
                                          str(phageInfo_handle.phageName) + ';' + str(spacer_id) + ';' + \
                                          str(missmatch) + ';' + str(start) + '-' + str(end) + '-' + '(%s);' % (order) + \
                                          str(flag) + '\n')
            interactionsfile_handle.flush()
            log_file.write('%d/%d  asn finished!\n' %(index,filenumber)  )
            log_file.flush()



if __name__ == "__main__":
    save_path = '/zrom/tmp/sql/'
    asnDir = '/zrom/tmp/archaeadbinteraction/'
    now_time = time.strftime('%Y-%m-%d', time.localtime(time.time()))
    phage_file_path = '/zrom/jobs/viral/'
    idToNameDic = get_id_to_name('/zrom/jobs/list/bacterialId')
    idToGenus_specisDic = get_genus_specis('/zrom/jobs/list/idToSpecisGenus.txt')
    interactionsfile_handle = open(save_path + "interation" + now_time + '.csv', 'w')
    hitproteinfile_handle = open(save_path + 'hitproteins' + '.csv', 'w')
    seqfile_handle = open(save_path + "sequences" + '.csv', 'w')
    interactionsfile_handle.write('bacteria_id;bacterial_name;specie;genus;phage_id;phage_name;spacer_id;missmatch;phage_position;flag\n')
    seqfile_handle.write('bac_id;phage_id;spacer_id;spacer_sequence;hit_sequence;long_hit_sequence;coverage\n')
    hitproteinfile_handle.write('bac_id;phage_id;spacer_id;protein_id;protein_position;protein_product;protein_sequence')
    log_file = open('/zrom/tmp/log/parserasn.log','w')
    paserAsn(asnDir)



