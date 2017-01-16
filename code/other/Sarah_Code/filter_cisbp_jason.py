import sys as sys

def filter_my_list(ffile,infile,outfile):
    o1 = open(ffile)
    o2 = open(infile)
    o3 = open(outfile,'w')
    lines = o1.readlines()
    motifs = {}
    for l in lines:
        s = l.strip().split('\t')
        motifs[s[0]] = s[1]

    for l in o2:
        s = l.strip().split('\t')
        if s[1] in motifs:
            o3.write(l)

    o3.close()
    o2.close()
    o1.close()

if __name__=="__main__":
    #sys.argv[1] has the name of the file which contains “good” motifs
    #    should have 1 motif ID per line
    # sys.argv[2] has the name of the file that is to be filtered for file in
    #     sys.argv[1]
    # sys.argv[3] is the name of the output file
    filter_my_list(sys.argv[1],sys.argv[2],sys.argv[3])



