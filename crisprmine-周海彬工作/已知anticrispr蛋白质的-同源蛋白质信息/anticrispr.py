#/usr/bin/python 
#!-*- coding:utf-8 -*-
import re
from Bio.Blast import NCBIXML
from Bio import Entrez
from Bio import SeqIO
import os
def mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)
def paser(filename):
    result_handle=open(filename)
    pattern=re.compile(r'\[(.*?)\]',re.S|re.I)
    blast_records=NCBIXML.parse(result_handle)
    for record in blast_records:
        queryName=record.query
        for index,al in enumerate(record.alignments):
            subjectNames=al.title
            for hsp in al.hsps:
                evalue=hsp.expect
                count=0
                for a in hsp.match:
                    if a=='+' or a==' ':
                        count+=1
                identy=float(hsp.positives-count)/(len(hsp.query))
                subjectQuery=hsp.sbjct
                if evalue<0.01 and identy>0.3 and al.qcovs>80:
                    queryId=queryName.split(' ')[0]
                    subjectNames=subjectNames.split(' >')
                    for subjectName in subjectNames:
                        elements=subjectName.split('|')
                        gi=elements[1]
                        geneType=elements[2]
                        proteinId=elements[3]
                        match=pattern.search(elements[4])
                        if match:
                            matchName=match.group()
                        else:
                            matchName='unkown'
                        if geneType=='gb':
                            a=HomologousProtein(anticrispr=queryId,proteinId=proteinId,geneIdType=geneType)
                            a.getGeneId()
                            try:
                                if matchName.find('phage')!=-1:
                                    a.getProteinsFromGene(savePath=phageSavePath,genePath=phageGenePath)
                                else:
                                    a.getProteinsFromGene(savePath=bacSavePath,genePath=bacGenePath)
                            except:
                                continue
                            print('>')
                            print(queryName)
                            print('|||||||||')
                            print(proteinId+'-'+geneType+'-'+a.geneId+'-'+matchName)
class HomologousProtein(object):
    def __init__(self,anticrispr,proteinId,geneIdType):
        Entrez.email='zhou13146@gmail.com'
        self.anticrispr=anticrispr
        self.proteinId=proteinId
        self.geneId=None
        self.geneIdType=geneIdType
    def getGeneId(self):
        if self.geneIdType=='gb':
            handle=Entrez.efetch(db='protein',id=self.proteinId,rettype='gb',retmode='text')
            assessionPattern=re.compile(r'DBSOURCE(.*?)KEYWORDS',re.S|re.I)
            match=assessionPattern.search(handle.read())
            if match:
                geneId=match.group().split('\n')[0].split(' ')[-1]
                self.geneId=geneId
            else:
                print(self.proteinId+'none find~')
    def getProteinsFromGene(self,savePath,genePath=None):
        if self.geneId!=None:
            savefile=open(savePath+self.geneId+'.faa','w')
            handle=Entrez.efetch(db="nucleotide",id=self.geneId,rettype="gb",retmod="text")
            if genePath!=None:
                SeqIO.convert(handle,'genbank',genePath+self.geneId+'.fna','fasta')
                handle.close()
                handle=Entrez.efetch(db="nucleotide",id=self.geneId,rettype="gb",retmod="text")
            records=list(SeqIO.parse(handle,'gb'))
            handle.close()
            for record in records:
                for feature in record.features:
                    if feature.type=='CDS':
                        location=feature.location
                        if str(location).find('+')!=-1:
                            direction='+'
                        elif str(location).find('-')!=-1:
                            direction='-1'
                        product=feature.qualifiers['product'][0] if feature.qualifiers.has_key('product') else 'unkown'
                        proteinId=feature.qualifiers['protein_id'][0] if feature.qualifiers.has_key('protein_id') else 'unkown'
                        translation=feature.qualifiers['translation'][0] if feature.qualifiers.has_key('translation') else 'unkown'
                        savefile.write('>ref'+'|'+str(proteinId)+'|'+str(location)+'|'+str(product)+'\n')
                        if translation[-1]=='\n':
                            savefile.write(translation)
                        else:
                            savefile.write(translation+'\n')
            savefile.close()
if __name__=='__main__':
    bacSavePath='bacfaa/'
    bacGenePath='bacfna/'
    phageSavePath='phagefaa/'
    phageGenePath='phagefna/'
    mkdir(bacGenePath)
    mkdir(bacSavePath)
    mkdir(phageGenePath)
    mkdir(phageSavePath)
    paser('anticrispr.xml')
    '''
    a=HomologousProtein(anticrispr='YP_007392342.1',proteinId='EZO98398.1',geneIdType='gb')
    a.getGeneId()
    a.getProteinsFromGene(savePath='',genePath=' ')
    '''
