# this script gets the fasta sequences for all peaks in the input file

import os as os
import sys as sys


def write_fasta(ifile,ofile):

    chrD = {}
    for i in range(0,21):
        print i
        o = open('/Users/SaraM/Dropbox/Collaborations/ImmGenATAC/SMdata/mm10/Mus_musculus.GRCm38.dna.chromosome.' + str(i+1) + '.fa')
        l = o.readline()
        
        s = ''
        for j in o:
            s += j.strip()
        chrD['chr'+str(i+1)] = s
        o.close()       

    o = open(ifile)
    t1 = open(ofile+'.chr1-5.txt','w')
    t2 = open(ofile+'.chr6-12.txt','w')
    t3 = open(ofile+'.chr13-19.txt','w')
    t4 = open(ofile+'.chrX-Y.txt','w')

    for l in o:
        s = l.strip().split('\t')

        if s[0]=='chrX':
            s[0] = 'chr20'
        if s[0]=='chrY':
            s[0] = 'chr21'

        try:
            mystring = chrD[s[0]][int(s[1]):int(s[2])]
        except:
            continue

        if int(s[0][3:])<6:
            o2 = t1
        elif int(s[0][3:])>5 and int(s[0][3:])<13:
            o2 = t2
        elif int(s[0][3:])>12 and int(s[0][3:])<20:
            o2 = t3
        else:
            o2 = t4
        o2.write('>'+s[0]+'-'+s[1]+':'+s[2]+':'+s[3]+'\n')
        for i in range(0,len(mystring)):
            if (i) % 50 == 0 and i>0:
                o2.write('\n')
            o2.write(mystring[i])
        o2.write('\n')
    o2.close()
    
        

if __name__=="__main__":
    write_fasta(sys.argv[1],sys.argv[2])
