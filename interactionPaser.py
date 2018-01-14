#! /usr/bin/env python
from Bio import Seq
import  os
from Bio import SeqIO
from Bio import Entrez
import threading
from time import sleep
from xml.etree import ElementTree
Entrez.email = "coldsummerday1211@gmail.com"
def get_file_path_hash(refseqid):
    number = 0
    for a in refseqid:
        number+=ord(a)
    saveNumber = number%150
    return str(saveNumber)

def get_spc_sequence(spcfile,spacerid):
    spacerid=">" + spacerid
    with open(spcfile,'r') as  file_handle:
        lines = file_handle.readline()
    for index,line in enumerate(lines):
        line = line[0:len(line)-1]
        if line == spacerid:
            seq = lines[index+1]
            if seq[-1] == '\n':
                seq = seq[0:len(seq)-1]
        return seq


class BacterialAsmParser(object):
    def __init__(self,bacterialId,asmSavePath):
        self.bacterialPath ='/zrom/jobs/bacteria/'
        self.bacterialId =bacterialId
        self.asmPath = asmSavePath
        self.interaction_list=[]
        self.empty_flag = False
        self.parserFile()
        self.seq = ''

        self.getBacterialSeq()

    def init_all(self):
        a=1

    def get_species_genus(self):
        if not  os.path.exists(self.bacterialPath + \
                  get_file_path_hash(self.bacterialId) + '/' + \
                  self.bacterialId + '/' \
                  + self.bacterialId+'.gb'):
            socketBuffer=''
            while socketBuffer=='':
                sleep(0.15)
                handle = Entrez.efetch(db="nucleotide", id=self.bacterialId, rettype="gb", retmode="text")
                socketBuffer = handle.read()
            with open(self.bacterialPath + \
                  get_file_path_hash(self.bacterialId) + '/' + \
                  self.bacterialId + '/' \
                  + self.bacterialId+'.gb','w') as file_handle:
                file_handle.write(socketBuffer)
        records = list(SeqIO.parse(self.bacterialPath + \
                  get_file_path_hash(self.bacterialId) + '/' + \
                  self.bacterialId + '/' \
                  + self.bacterialId+'.gb', "genbank"))
        taxonomylist =records[0].annotations['taxonomy']



    def getBacterialSeq(self):
        with open(self.bacterialPath + \
                  get_file_path_hash(self.bacterialId) + '/' + \
                  self.bacterialId + '/' \
                  + self.bacterialId+'.spc','r') as file_handle:
            lines = file_handle.readlines()
        allline = ''
        for line in lines:
            if line[0]=='>':
                continue
            line=line[0:len(line)-1]
            allline+=line
        self.seq = allline

    def parserFile(self):
        with open(self.asmPath+self.bacterialId+'.asn') as file_handle:
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
            end = int (elements[9])
            length = alength / identy
            if length < 22:
                self.interaction_list.append \
                    ({'spacerid':spacerid,'phageid':phageid,\
                      'missmatch':missmatch,'start':start,'end':end})
                continue
            if missmatch < ((length - 22) ** 0.5):
                self.interaction_list.append \
                    ({'spacerid': spacerid, 'phageid': phageid,\
                      'missmatch': missmatch,'start':start,'end':end})



class phageInfo(object):
    def __init__(self,phageID,phagefilePath):
        self.phageID = phageID
        self.phagefilePath = phagefilePath
        self.proteinList = []
        self.phageName = None
        self.phageSeq = ''
        self.get_phage_seq()
        self.parserGenbankFile()


    def get_phage_seq(self):
        with open(self.phagefilePath+self.phageID+'.fna','r') as file_handle:
            lines = file_handle.readlines()
        alllines=''
        for line in lines:
            if line[0]=='>':
                continue
            line=line[0:len(line)-1]
            alllines+=line
        self.phageSeq = alllines

    def get_position_seq(self,start,end):
        return self.phageSeq[start:end]


    def parserGenbankFile(self):
        if not  os.path.exists(self.phagefilePath+self.phageID+'.gb'):
            try:
                sleep(0.05)
                handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=self.phageID)
                with open(self.phagefilePath+self.phageID+'.gb', 'w') as file_handle:
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
                    protein = phageProtein(proteinId,location,product,translation)
                    protein.parserLocation()
                    self.proteinList.append(protein)

    def get_interaction_protein(self,start,end,order):
        for protein in self.proteinList:
            if order=='+' and protein.order == '+':
                if start > protein.start and end < protein.end:
                    return protein
            elif order == '-' and protein.order == '-':
                if start > protein.start and end < protein.end:
                    return protein
        return None



class phageProtein(object):
    def __init__(self,proteinID,location,product,seq):
        self.proteinID = proteinID
        self.location = location
        self.product = product
        self.seq = seq

    def parserLocation(self):
        start = int(self.location.split(':')[0].split('[')[1])
        end =   int(self.location.split(':')[1].split(']')[0])
        self.order = self.location.split('(')[1].split(')')[0]
        if start > end:
            temp = start
            start = end
            end = temp
        self.start = start
        self.end = end




if __name__=="__main__":
    print(get_file_path_hash("NC_023738"))
