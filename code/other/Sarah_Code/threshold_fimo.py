import sys as sys

def mythreshold(infile,thr,outfile):

    o = open(infile)
    o2 = open(outfile,'w')
    for l in o:
        s = l.strip().split('\t')
        if float(s[3])<float(thr):
            o2.write(l)


    o2.close()
    o.close()

if __name__=="__main__":
    mythreshold(sys.argv[1],sys.argv[2],sys.argv[3])


