#!/usr/bin/python
#!-*-coding:utf-8-*-
from Bio.Blast.Applications import  NcbiblastnCommandline
import subprocess
import os
from os.path import getsize
from subprocess import call
import argparse
class Interation(object):
    def __init__(self,bacteriaId,spcFileName):
        self.bacteriaId=bacteriaId
        self.spcFileName=spcFileName
        self.db='/zrom/jobs/phagefna/phagefna'
        self.savepath='/zrom/tmp/tmpasn/'
    def getInteraction(self):
        blast_cline="blastn -query %s -db %s -outfmt 6 -evalue 0.01 -word_size 18 -out %s" % (self.spcFileName,self.db,self.savepath+self.bacteriaId+'.asn') 
        blast_cline=NcbiblastnCommandline(query=self.spcFileName,db=self.db,evalue=1e-2,outfmt=6,word_size=18,out=self.savepath+self.bacteriaId+'.asn')
        try:
            blast_cline()
        except Exception as e:
            pass 
        result_file=open(self.savepath+self.bacteriaId+'.res','w')
        if os.path.exists(self.savepath+self.bacteriaId+'.asn'):
            size=getsize(self.savepath+self.bacteriaId+'.asn')
            if size==0:
                result_file.write("None")
                result_file.close()
                call("rm "+self.savepath+self.bacteriaId+".asn",shell=True)
                return
            else:
                
                with open(self.savepath+self.bacteriaId+'.asn') as file_handle:
                    lines=file_handle.readlines()
                interactions=[]
                for line in lines:
                    elements=line.split('\t')
                    spacerId=elements[0]
                    subject=elements[1]
                    identy=float(elements[2])
                    alength=float(elements[3])
                    missmatch=float(elements[4])
                    position=elements[8]+'-'+elements[9]
                    info={'subject':subject,'spacerId':spacerId,'identity':identy,'len':int(alength),'position':position}
                    if alength<22:
                        if info not in interactions:
                            interactions.append(info)        
                        continue
                    if missmatch< ((alength-22) ** 0.5):
                        if info not in interactions:
                            interactions.append(info)    
            call("rm "+self.savepath+self.bacteriaId+".asn",shell=True)
            for a in interactions:
                seq=getsequen(self.spcFileName,a['spacerId'])
                cov=alength/len(seq)
                result_file.writelines("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(a['subject'],a['spacerId'],str(a['identity']),cov,a['len'],a['position'],seq))
            result_file.close()
            return
def getsequen(spcfile,spacerId):
    spacerid=">"+spacerId
    with open(spcfile,'r') as file_handle:
        lines=file_handle.readlines()
    for index,line in enumerate(lines):
        line=line[0:len(line)-1]
        if line==spacerid:
            seq=lines[index+1]
            if seq[-1]=='\n':
                seq=seq[0:len(seq)-1]
            return seq
def getlist(filename,num):
    with open(filename,'r') as file_handle:
        lines=file_handle.readlines()
    returnlist=[]
    for line in lines:
        line=line[0:len(line)-num]
        returnlist.append(line)
    return returnlist
def init():
    parser=argparse.ArgumentParser()
    parser.add_argument('-s','--spc',type=str,help="spcfile include path")
    args=parser.parse_args()
    if args.spc:
        return args.spc
    else:
        print('input the true file')
        exit(1)
if __name__=="__main__":
    spcfile=init()
    geneId=spcfile.split('/')[-1]
    geneId=geneId[0:len(geneId)-4]
    inter=Interation(geneId,spcfile)
    inter.getInteraction()
    '''
    lists=getlist('spclist',5)
    for spc in lists:
        a=Interation(spc,'/home/admin01/spc/%s.spc' %(spc))
        a.getInteraction()
    '''
                        

        
            
            
        