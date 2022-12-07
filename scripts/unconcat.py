import sys
import glob
import os

#this scripts is unconcatenating Hit-scIsoseq reads (after performing a blast search)

def read_blast(path):

    blast={}
    with open(path, 'r') as m7:
        lines = m7.readlines()
        for i in range(len(lines)):
            if lines[i][0]=='#':continue
            [QNAME,SUBJECT,_,AL_LENGTH,_,_,START,END,_,_,EVAL,_] = lines[i].split()
            if float(EVAL)<0.005:
                try:
                    blast[QNAME].append([int(START),int(END)])
                except KeyError:
                    blast[QNAME]=[[int(START),int(END)],]

    for key in blast:
        blast[key]=sorted(blast[key])
    
    return blast


def main(samfile,m7file,sample):
    
    blast = read_blast(m7file)
    junk = []
    is0 = {k: [] for k in blast.keys()}
    p=1
    
    with open(samfile, 'r') as sam, open(sample, 'w') as newsam:

        lines = sam.readlines()

        for line in lines:
            
            if line[0]=='@':
                newsam.write(line)
                continue
            
            [QNAME,FLAG,RNAME,POS,MAPQ,CIGAR,_,_,_,SEQ,QUAL] = line.split()[:11]

            MAX = len(SEQ)
            newline=line.split()
            ZMW=newline[-1]
            
            try:
                positions = blast[QNAME]
            except KeyError:
                newline[-1]=':'.join(ZMW.split(':')[:2]+[str(p)+'\n'])
                newsam.write("\t".join(newline))
                junk.append(QNAME)
                p+=1
                continue
                
            try:
                osef=positions[1]
            except IndexError:
                newline[-1]=':'.join(ZMW.split(':')[:2]+[str(p)+'\n'])
                newsam.write("\t".join(newline))
                junk.append('zmw:i:'+str(p))
                p+=1
                continue

            if positions[0][1] >= 400:
                junk.append(QNAME)

            if positions[-1][1] <= MAX-400:
                junk.append(QNAME)
    
            
            i=0
            while positions[i][1]<=MAX:
                try:
                    if positions[i+1][1] - positions[i][0] > 400:
                        is0[QNAME].append((positions[i][0],positions[i+1][1]))
                        i+=1
                    else:
                        i+=1
                except IndexError:
                    break

            is0[QNAME]=sorted(list(set(is0[QNAME])))
            
            for i in range(len(is0[QNAME])):
                (start,stop) = is0[QNAME][i]
                newline[9:11]=[SEQ[start-1:stop], QUAL[start-1:stop]]
                newline[-1]=':'.join(ZMW.split(':')[:2]+[str(p)+'\n'])
                newsam.write("\t".join(newline))
                p+=1
                
    print("There were %i reads containing some trash..." % len(junk))


# calling main function with the command line (input = path to sam) 
if __name__ == '__main__': 
    if len(sys.argv) != 4:
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
