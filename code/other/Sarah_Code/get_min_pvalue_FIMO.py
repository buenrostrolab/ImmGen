import sys as sys
import os as os

def write_max_score(ifile,ofile):
    o1 = open(ifile)
    o2 = open(ofile,'w')
    
    peak = {}
    for l in o1:
        s = l.strip().split('\t')
        if (s[0],s[1]) in peak:
            if float(s[3])<peak[s[1]][1]:

                peak[(s[1],s[0])] = (float(s[2]),float(s[3]))
        else:
            peak[(s[1],s[0])] = (float(s[2]),float(s[3]))
    
    
    for k in peak.keys():
        o2.write(k[0]+'\t'+k[1]+'\t'+str(peak[k][0])+'\t'+str(peak[k][1])+'\n')
    o2.close()

if __name__=="__main__":
    # sys.argv[1] contains the name of the FIMO output file, which contains
    # p-values as 4th column and peak and motif IDs as first and second columns
    # sys.argv[2] is the name of the output file
    write_max_score(sys.argv[1],sys.argv[2])

     
