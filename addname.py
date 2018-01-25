


def get_file_path_hash(refseqid):
    number = 0
    for a in refseqid:
        number += ord(a)
    saveNumber = number % 150
    return str(saveNumber)


def getqfastalist(file):
    with open(file,'r') as file_handle:
        lines=file_handle.readlines()
    list=[]
    for line in lines:
        line = line.split('\t')[0]
        list.append(line)
    return list

def getname(id):
    with open('/zrom/jobs/bacteria/'+get_file_path_hash(id)+'/%s/%s.fna' %(id,id) \
            ,'w') as file_handle:
        line = file_handle.readline()
    name = line.split('|')[-1]
    return name


if __name__ == '__main__':
    newlist = getqfastalist('/zrom/jobs/list/cripsrarchaea')

    savefile = open('/zrom/jobs/list/bacterialId','a')
    for id in newlist[:50]:
        name = getname(id)
        print(id,name)
        #savefile.write(id+'---- '+name+'\n')
    savefile.close()
