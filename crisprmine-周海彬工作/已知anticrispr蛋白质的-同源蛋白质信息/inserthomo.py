#!/usr/bin/python 
#-*-coding:utf-8-*-
import MySQLdb
def parser(filename,con,cur):
    with open(filename,'r') as file_handle:
        lines=file_handle.readlines()
    allline=''
    for line in lines:
        allline+=line
    elements=allline.split('>')
    elements.remove("")
    for index,hit in enumerate(elements):
        antiline,homoline=hit.split("|||||||||")
        antiline=antiline[1:len(antiline)-1]
        homoline=homoline[1:len(homoline)-1]
        antiproteinId=antiline.split(' ')[0]
        homoelemnts=homoline.split('-')
        proteinId=homoelemnts[0]
        genetype=homoelemnts[1]
        geneid=homoelemnts[2]
        genename=homoelemnts[3]
        
        sql="INSERT INTO homologyprotein(proteinId, genetype, geneid, genename, antiproteinId)VALUES('%s', '%s', '%s', '%s', '%s');" %(proteinId,genetype,geneid,genename,antiproteinId)
        print(sql)
        
    cur.close()
    con.close()
def connect():
    con=MySQLdb.connect("localhost","root","zhouhaibin1211","bio")
    cur=con.cursor()
    
    return con,cur
        
if __name__=="__main__":
    con,cur=connect()
    parser('antiHomo.txt',con,cur)