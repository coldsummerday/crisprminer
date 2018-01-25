#!/usr/bin/python
#!-*-coding:utf-8:-*-
import re



def parser(geneid):
    global cddfile
    global proteinsfile
    exitPattern=re.compile(r'-----',re.I|re.S)
    with open('proteins/%s.ph' % geneid,'r') as file_handle:
        lines=file_handle.readlines()
    allline=''
    for line in lines:
        allline+=line
    proteins=allline.split('>>>')
    for protein in proteins:
        math=exitPattern.search(protein)
        if math!=None:
            if '' in  proteins:
                proteins.remove('')
            hits=protein.split("-----")
            lines=hits[0]
            hits.remove(lines)
            for a in hits:
                if a=='\n':
                    hits.remove(a)
            lines=protein.split('\n')
            if '' in lines:
                lines.remove('')
            elements=lines[0].split('|')
            query=lines[1]

            proteinId=elements[1]
            product=elements[3]
            position=elements[2]
            proteinsql="INSERT INTO proteins(proteinId, geneid, product, geneposition, query)VALUES('%s', '%s', '%s', '%s', '%s');" %(proteinId,geneid,product,position,query)
            proteinsfile.write(proteinsql+'\n\n')
            for hit in hits:
                if hit=="\n" or hit=="":
                    continue
                templines=hit.split("\n")
                templine=templines[1]
                elements=templine.split("	 ")
                evalue=elements[0].split(":")[1]
                identy=elements[1].split(":")[1]
                cddline=templines[2]
                index=cddline.find(".")
                description=cddline[index+1:]
                cddnumline=cddline[0:index]
                shortname=cddnumline.split(",")[1]
                if len(cddnumline.split(","))>2:
                    description=cddnumline.split(",")[2]+'.'+description
                description=description.replace("'",'"')
                cddid=cddnumline.split(",")[0].split(" ")[-1]
                cddsql="INSERT INTO proteindomain(proteinId,geneid, proteindomain, cdd, cdddescription, evalue, identy)VALUES('%s','%s' ,'%s', '%s', '%s', '%s', '%s');" %(proteinId,geneid,shortname,cddid,description,evalue,identy)
                cddfile.write(cddsql+'\n\n')
        else:
            if protein=="":
                continue
            lines=protein.split('\n')
            if '' in lines:
                lines.remove('')
            elements=lines[0].split('|')
            query=lines[1]
            proteinId=elements[1]
            product=elements[3]
            position=elements[2]
            proteinsql="INSERT INTO proteins(proteinId, geneid, product, geneposition, query)VALUES('%s', '%s', '%s', '%s', '%s');" %(proteinId,geneid,product,position,query)
            cddsql="INSERT INTO proteindomain(proteinId, geneid , proteindomain, cdd, cdddescription, evalue, identy)VALUES('%s','%s' ,'%s', '%s', '%s', '%s', '%s');" %(proteinId,geneid,'NULL','NULL','NULL','NULL','NULL')
            proteinsfile.write(proteinsql+'\n\n')
            cddfile.write(cddsql+'\n\n')
def getlist(filename,number):
    with open(filename,'r') as file_handle:
        lines=file_handle.readlines()
    list=[]
    for line in lines:
        line=line[0:len(line)-number]
        list.append(line)
    return list


if __name__=="__main__":
    phlist=getlist("list/phlist",4)
    proteinsfile=open("insertproteins.sql",'w')
    cddfile=open("insertproteindomain.sql",'w')
    for ph in phlist:
        parser(ph)
    proteinsfile.close()
    cddfile.close()
    
        
    