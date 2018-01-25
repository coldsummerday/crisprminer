#!/usr/bin/python

# -*- coding:utf-8 -*-


from Bio.Blast.Applications import NcbiblastpCommandline
import subprocess
import os
from os.path import getsize
from subprocess import call
import argparse
from Bio import SeqIO
from Bio.Blast import NCBIXML
class anticirspr(object):
    def __init__(self,proteinfile,geneid,savePath):
        self.proteinfile=proteinfile
        self.savePath=savePath
        self.geneid=geneid
        self.dbpath='/zrom/jobs/anticrisprprotein/anticrisprprotein'
    def blastp(self):
        blast_cline=NcbiblastpCommandline(query=self.proteinfile,
                                        db=self.dbpath,out=self.savePath+self.geneid+'.xml',
                                        qcov_hsp_perc=80,
                                        outfmt=5,
                                        evalue=1e-2)
        a=''
        try:
	    #print(blast_cline)
            blast_cline()
        
        except Exception as es:
            a='error'
        
        result_file=open(self.savePath+self.geneid+'.anti','w')
        if os.path.exists(self.savePath+self.geneid+'.xml'):
            flag=True
        else:
            flag=False
        hitlist=[]
        if flag:
            blast_records=NCBIXML.parse(open(self.savePath+self.geneid+'.xml'))
            for record in blast_records:
                queryName=record.query
                for index,al in enumerate(record.alignments):
                    subjectName=al.title
                    for hsp in al.hsps:
                        seq=getseqs(self.proteinfile,queryName)
                        qcovs=len(hsp.sbjct)/len(seq)
                        evalue=hsp.expect
                        count=0
                        for a in hsp.match:
                            if a=='+' or a==' ':
                                count+=1
                        identy=abs(float(hsp.positives-count)/(len(hsp.query)))
                        subjectQuery=hsp.sbjct
                    if evalue<0.01 and identy>0.3 and qcovs>0.8:
                        hitlist.append({"protein":queryName,"hitanti":subjectName,"evalue":str(evalue),'identy':str(identy),'cov':str(qcovs),'proteinquery':seq})
        if len(hitlist)>0:
            for line in hitlist:
                result_file.write(line['protein']+'\t'+line['hitanti']+'\t'+line['evalue']+'\t'+line['identy']+'\t'+line['cov']+'\t'+line['proteinquery']+'\n')
        else:
            result_file.write('None')
    

def getseqs(filename,proteinid):
    spacerid=">"+proteinid
    with open(filename,'r') as file_handle:
        lines=file_handle.readlines()
    for index,line in enumerate(lines):
        line=line[0:len(line)-1]
        if line==spacerid:
            seq=lines[index+1]
            if seq[-1]=='\n':
                seq=seq[0:len(seq)-1]
            return seq
def init():
    parser=argparse.ArgumentParser()
    parser.add_argument('-f','--file',type=str,help="proteinfile include path")
    args=parser.parse_args()
    if args.file:
        return args.file
    else:
        print('input the true file')
        exit(1)
if __name__=="__main__":
    filename=init()
    geneid=filename.split('/')[-1]
    geneid=geneid[0:len(geneid)-4]
    a=anticirspr(filename,geneid,'/zrom/tmp/tmpanti/')
    a.blastp()
        
        

        
