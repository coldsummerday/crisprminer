#!/usr/bin/env python
# -*- coding:utf-8 -*- 
import  os
import subprocess
from os.path import getsize
from Bio.Blast.Applications import NcbiblastnCommandline

def gethash(refseqid):
    number=0
    for a in refseqid:
        number+=ord(a)
    saveNumber=number%150
    return str(saveNumber)
def makedir(Path):
    if not os.path.exists(Path):
        os.mkdir(Path)
    #该函数作用为删除空文件
def rmfile(folder,refseqaccession):
    try:
    	size=getsize('%s/%s' %(folder,refseqaccession))
    	if size==0:
            subprocess.call(" rm %s/%s" % (folder,refseqaccession),shell=True)
    except:
        pass
def tryblast(refseqaccession):
    global BASE_DIR
    global savepath
    blast_cline=NcbiblastnCommandline(query=BASE_DIR+'jobs/bacteria/%s/%s/' % (gethash(refseqaccession),refseqaccession)+refseqaccession+".spc",
                                     db='/zrom/jobs/phagefna/phagefna',
                                      word_size=20,out=savepath+refseqaccession+".asn",
                                      outfmt="6",evalue=1e-2,
                                      )
    try:
        #进行Blast，并保存日志信息
        blast_cline()
        log_handle.write(str(count)+'/'+str(numbers)+'\t'+refseqaccession+'\t finished\n')
        log_handle.flush()
        
    except:
        log_handle.write(str(count)+'/'+str(numbers)+'\t'+refseqaccession+'\terror\n')
        log_handle.flush()
        error_handle.write(refseqaccession+'\n')
        error_handle.flush()
def getlist(filename,number):
    with open(filename,'r') as file_handle:
        lines=file_handle.readlines()
    list=[]
    for line in lines:
        line=line[0:len(line)-number]
        list.append(line)
    return list
def getqfastalist(file):
    with open(file,'r') as file_handle:
        lines=file_handle.readlines()
    list=[]
    for line in lines:
        line = line.split('\t')[0]
        list.append(line)
    return list
if __name__=="__main__":
    BASE_DIR='/zrom/'
    savepath=BASE_DIR+'tmp/archaeadbinteraction/'
    log_handle=open(BASE_DIR+'tmp/log/archaeadbblastinteraction','w')
    makedir(savepath)
    qfastalist=getqfastalist('/zrom/jobs/list/cripsrbacteria')
    numbers=len(qfastalist)
    count=0
    error_handle=open('/zrom/tmp/error/archaeadbblastinteraction','w')
    for qfasta in qfastalist:
        count+=1
        tryblast(qfasta)
       # rmfile(BASE_DIR+'tmp/dbinteraction',qfasta+'.asn')  
    error_handle.close()
    log_handle.close()
                
            
            
   


#!/usr/bin/env python
# -*- coding:utf-8 -*- 
import  os
import subprocess
from os.path import getsize
from Bio.Blast.Applications import NcbiblastnCommandline

def gethash(refseqid):
    number=0
    for a in refseqid:
        number+=ord(a)
    saveNumber=number%150
    return str(saveNumber)
def makedir(Path):
    if not os.path.exists(Path):
        os.mkdir(Path)
    #该函数作用为删除空文件
def rmfile(folder,refseqaccession):
    try:
    	size=getsize('%s/%s' %(folder,refseqaccession))
    	if size==0:
            subprocess.call(" rm %s/%s" % (folder,refseqaccession),shell=True)
    except:
        pass
def tryblast(refseqaccession):
    global BASE_DIR
    global savepath
    blast_cline=NcbiblastnCommandline(query=BASE_DIR+'jobs/bacteria/%s/%s/' % (gethash(refseqaccession),refseqaccession)+refseqaccession+".spc",
                                     db='/zrom/jobs/phagefna/phagefna',
                                      word_size=20,out=savepath+refseqaccession+".asn",
                                      outfmt="6",evalue=1e-2,
                                      )
    try:
        #进行Blast，并保存日志信息
        blast_cline()
        log_handle.write(str(count)+'/'+str(numbers)+'\t'+refseqaccession+'\t finished\n')
        log_handle.flush()
        
    except:
        log_handle.write(str(count)+'/'+str(numbers)+'\t'+refseqaccession+'\terror\n')
        log_handle.flush()
        error_handle.write(refseqaccession+'\n')
        error_handle.flush()
def getlist(filename,number):
    with open(filename,'r') as file_handle:
        lines=file_handle.readlines()
    list=[]
    for line in lines:
        line=line[0:len(line)-number]
        list.append(line)
    return list
def getqfastalist(file):
    with open(file,'r') as file_handle:
        lines=file_handle.readlines()
    list=[]
    for line in lines:
        line = line.split('\t')[0]
        list.append(line)
    return list
if __name__=="__main__":
    BASE_DIR='/zrom/'
    savepath=BASE_DIR+'tmp/archaeadbinteraction/'
    log_handle=open(BASE_DIR+'tmp/log/archaeadbblastinteraction','w')
    makedir(savepath)
    qfastalist=getqfastalist('/zrom/jobs/list/cripsrarchaea')
    numbers=len(qfastalist)
    count=0
    error_handle=open('/zrom/tmp/error/archaeadbblastinteraction','w')
    for qfasta in qfastalist:
        count+=1
        tryblast(qfasta)
       # rmfile(BASE_DIR+'tmp/dbinteraction',qfasta+'.asn')  
    error_handle.close()
    log_handle.close()
                
            
            
   



