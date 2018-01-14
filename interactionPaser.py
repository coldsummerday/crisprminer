#! /usr/bin/env python

import  os

def get_file_path(refseqid):
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


if __name__=="__main__":
    print(get_file_path("NZ_NWLE01000020"))
