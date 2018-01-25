#!/usr/bin/python
#-*-coding:utf-8 -*-

import sys
from Bio.Blast.Applications import NcbirpsblastCommandline
import os
import subprocess
import threading
import argparse
def rpsblast(phageid):
    global fnaPath
    global savePath
    blast_cline="rpsblast -query %s/%s.faa -comp_based_stats 0 -evalue 0.01 -seg no -outfmt 5 -num_threads 4 -db /zrom/z-tools/data/cdd/Cdd -out %s.xml" %(fnaPath,phageid,savePath+phageid)
    '''
    blast_cline=NcbirpsblastCommandline(query=fnaPath+phageid+'/'+phageid+'.afaa',
                                        evalue=0.01,seg='no',outfmat=7,
                                        db="/zrom/z-tools/data/cdd/Cdd | rpsbproc",
                                        d='/zrom/z-tools/data/cdd',
                                       )
    '''
    #进行blast
    try:
        code=subprocess.call(blast_cline,shell=True)
        if code==0:
            return None
        else:
            return phageid
    except:
        return phageid
def getlist(filename,number):
    with open(filename,'r') as file_handle:
        lines=file_handle.readlines()
    list=[]
    for line in lines:
        line=line[0:len(line)-number]
        list.append(line)
    return list
def mkdir(directory):
    if  not os.path.exists(directory):
        os.mkdir(directory)
class tThread(threading.Thread):
    def __init__(self,taskl):
        threading.Thread.__init__(self)
        self.taskl=taskl
    def run(self):
        try:
            blasp_loop(self.taskl)
        except:
            print("thread end")
def blasp_loop(baclist):
    for bac in baclist:
        rpsblast(bac)
def initArgs():
    global threadnum
    paser=argparse.ArgumentParser()
    parser.add_argument('-n','--number',type=int,help='thread number',default=1)
    args=paser.parse_args()
    if args.number:
        threadnum=args.number
if __name__=="__main__":
    threadnum=1
    fnaPath='./bacfaa/'
    savePath='bacxmls/'
    mkdir(savePath)
    baclist=getlist('list/bacfaalist',5)
    taskmanager={}
    taskid=0
    for taskm in baclist:
        taskid+=1
        taskgid=taskid%threadnum
        if (taskgid not in taskmanager.keys()):
            taskmanager[taskgid]=[]
        taskmanager[taskgid].append(taskm)
    for taskg in taskmanager.keys():
        t1=tThread(taskmanager[taskg])
        t1.start()
    