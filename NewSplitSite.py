#!/usr/bin/pyton2.7

'''
    This program is used to extract the specific splitsite.
    Input: a junction file
    Output: five file(existed junction site and condition ABCD)

    Author: XiaolongZhang
    Date: 2017-03-20
'''

import re
import pdb
import sys
from operator import itemgetter

Chrom = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
         'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20',
         'chr21','chr22','chrX', 'chrY', 'chrM']


# Global Match for gene_id, transcipt_id and gene_name
geneID = re.compile('gene_id ".*?";')
transID = re.compile('transcript_id ".*?";')
gname = re.compile('gene_name ".*?";')


class Trans(object):
    def __init__(self, gname, tid, gid):
        self.gname = gname        # gene name for transcript
        self.tid   = tid          # transcript ID
        self.gid   = gid          # gene ID
        self.Elist = []           # Exon list for transcript

    def AddExon(self, exon):
        ''' Add exon tuple for Elist'''
        self.Elist.append(exon)


def usage(prog):
    print 'python %s <*.gtf> <*.junc>' %(prog[0])
    sys.exit(1)


def ReadJunc(juncfile, juncdict):
    ''' juncdict = {'chr1':[(16027,16607, 'rev', '7'), (),...], 
                    'chr2':[(),(), ...], ... }
    '''
    print 'Start to read the junction file ...'
    fp, prechr = open(juncfile, 'r'), 'dummy'

    for line in fp.readlines():
        llist = line.strip().split('\t')
        content = ( int(llist[1]), int(llist[2]), llist[3], llist[4] )

        if llist[0] in Chrom:
            if llist[0] == prechr:
                juncdict[prechr].append(content)
            else:
                prechr = llist[0]
                juncdict[prechr] = [content]
    fp.close()


def ReadGtf(gtffile, gtfdict):
    ''' gtfdict = {'chr1':[T1, T2, T3, ...], 'chr2':[T1,T2,T3, ...], ...}
        T1 = Trans(gname, tid, gid)
        T1.AddExon( (estar, eend), ... )
    '''
    print 'Start to read the gtf file ...'
    fp, prechr = open(gtffile, 'r'), 'dummy'
    firsttrans, T = True, Trans('dummy','dummy','dummy')

    for line in fp.readlines():
        if not line.startswith('#'):
            llist = line.strip().split('\t')
            if llist[0] == prechr:
                if llist[2] == 'transcript':
                    tid = transID.search(llist[8]).group()[15:-2]
                    gid = geneID.search(llist[8]).group()[9:-2]
                    genename = gname.search(llist[8]).group()[11:-2]
                    if firsttrans:
                        firsttrans =  False
                    else:
                        gtfdict[prechr].append( T )

                    T = Trans( genename, tid, gid )
                elif llist[2] == 'exon':
                    T.AddExon( (int(llist[3]), int(llist[4])) )
            else:
                if prechr != 'dummy':
                    gtfdict[prechr].append( T )
                prechr, gtfdict[prechr], firsttrans = llist[0], [], True

    gtfdict[prechr].append( T )
    fp.close()


def GetExonList(gtfdict, chrom, key=0):
    ''' key = 0: sorted by exon start [default]
        key = 1: sorted by exon end
    '''
    exonlist, glen = [], len(gtfdict[chrom])

    for i in xrange(glen):
        exon = gtfdict[chrom][i].Elist
        for j in xrange(len(exon)):
            exonlist.append(exon[j] + (i,j))

    if not key:
        return sorted(exonlist)
    else:
        return sorted(exonlist, key=itemgetter(1))


def GetJuncList(gtfdict, chrom):
    junclist, glen = [], len(gtfdict[chrom])

    for i in xrange(glen):
        exon = sorted(gtfdict[chrom][i].Elist)

        if len(exon) > 1:
            for j in xrange(len(exon)-1):
                junclist.append( (exon[j][1], exon[j+1][0]) )
    return sorted(junclist)
                


def BinSearch(pos, exonlist, locindex):
    low, high = 0, len(exonlist)-1 

    while (low <= high):
        mid = (low + high) / 2
        if exonlist[mid][locindex] < pos:
            low = mid + 1
        elif exonlist[mid][locindex] > pos:
            high = mid - 1
        else:
            return mid
    return (-1)


def LocatedSite(pos, exonlist, locindex=0):
    ''' locindex = 0: exonlist was sorted by start[default]
        locindex = 1: exonlist was soeted by end.
    '''
    targetexon = []

    posindex = BinSearch(pos, exonlist, locindex)
    if posindex == -1:
        return targetexon
    else:
        targetexon.append(exonlist[posindex])
        leftindex, rightindex = posindex -1, posindex +1

        # Search for the pos in the left site
        while (leftindex >= 0):
            if exonlist[leftindex][locindex] == pos:
                targetexon.append(exonlist[leftindex])
            else:
                break
            leftindex -= 1

        # Search for the pos in the right site
        maxindex = len(exonlist) - 1
        while (rightindex <= maxindex):
            if exonlist[rightindex][locindex] == pos:
                targetexon.append(exonlist[rightindex])
            else:
                break
            rightindex += 1
    return targetexon


def FiltJuncdic(gtfdict, juncdict):
    ''' junclist = [(1186,1442), (1186,1654), (1324,2341), ...]
    '''
    print 'Start to filt the existed junction ...'
    newjuncdict, existjunc = {}, []
    for chrom in juncdict.keys():
        newjuncdict[chrom] = []

    for c in Chrom:          # c: chromosome => chr1, chr2, ...
        if c in juncdict.keys():
            junclist = GetJuncList(gtfdict, c)
            for e in juncdict[c]:      # e: junctionline => (1,2,'rev','7'), ...
                temlist, status = LocatedSite(e[0], junclist), True
                if temlist:
                    for f in temlist:
                        if e[1] == f[1]:
                            status = False
                            existjunc.append( (c,) + e )
                            break
                    if status:
                        newjuncdict[c].append(e)
                else:
                    newjuncdict[c].append(e)

    with open("existjunc.txt", 'w') as efile:
        for line in existjunc:
            outlist = [ str(e) for e in line ]
            efile.write("%s\n" %('\t'.join(outlist)))

    return newjuncdict



def ConditionAC(gtfdictlist, chrom, juncline, e_exonlist, splitdict):
    ''' conditionA:              |         conditionC:
          ______                 |          _______________________
         |      |                |         |             |         |
    =====     ==*==     =====    |    =====     =====     =====     =====
         1   2     3   4     5   |         1   2     3   4     5   6
    '''
    loclist = LocatedSite(juncline[0], e_exonlist, 1)
    if loclist:
        locgtf = [ (gtfdictlist[e[2]], e[3]) for e in loclist ]
        for g in locgtf:
            # Condition A
            try:
                nextexon = (g[0].Elist[g[1]+1][0], g[0].Elist[g[1]+1][1])
                if juncline[1] > nextexon[0] and juncline[1] < nextexon[1]:
                    teminfo1 = chrom+'\t'+g[0].tid+'\t'+g[0].gname+'\t'+str(nextexon[0])+':'+str(nextexon[1])
                    teminfo2 = ':'.join([str(e) for e in juncline])
                    splitdict['A'].append(teminfo1+'\t'+teminfo2)
            except IndexError:
                pass

            # Condition C
            for i in range(g[1]+2, len(g[0].Elist)):
                texon = (g[0].Elist[i][0], g[0].Elist[i][1])
                if juncline[1] == texon[0]:
                    teminfo1 = chrom+'\t'+g[0].tid+'\t'+g[0].gname+'\t'+str(texon[0])+':'+str(texon[1])
                    teminfo2 = ':'.join([str(e) for e in juncline])
                    splitdict['C'].append(teminfo1+'\t'+teminfo2)
                
        

def ConditionBD(gtfdictlist, chrom, juncpair, s_exonlist, splitdict):
    ''' conditionB:              |         conditionD:
       ______                    |          
      |      |                   |         
    ==*==     =====     =====    |    =====     =====*****=====     =====
         1   2     3   4     5   |         1   2     3   4     5   6
    '''
    loclist = LocatedSite(juncpair[0][1], s_exonlist)
    if loclist:
        locgtf = [ (gtfdictlist[e[2]], e[3]) for e in loclist ]
        for g in locgtf:
            # Condition B
            try:
                preexon = (g[0].Elist[g[1]+1][0], g[0].Elist[g[1]+1][1])
                if juncpair[0][0] > preexon[0] and juncpair[0][0] < preexon[1]:
                    teminfo1 = chrom+'\t'+g[0].tid+'\t'+g[0].gname+'\t'+str(preexon[0])+':'+str(preexon[1])
                    teminfo2 = ':'.join([str(e) for e in juncpair[0]])
                    splitdict['B'].append(teminfo1+'\t'+teminfo2)
            except IndexError:
                pass

            # Condition D
            try:
                currexon = (g[0].Elist[g[1]][0], g[0].Elist[g[1]][1])
                nextexon = (g[0].Elist[g[1]-1][0], g[0].Elist[g[1]-1][1])
            except IndexError:
                continue

            if juncpair[1][0] == nextexon[1]:
                teminfo1 = chrom+'\t'+g[0].tid+'\t'+g[0].gname
                teminfo2 = str(currexon[0])+':'+str(currexon[1])+'/'+str(nextexon[0])+':'+ str(nextexon[1])
                teminfo3 = ':'.join([str(e) for e in juncpair[0]])+'/'+':'.join([str(e) for e in juncpair[1]])
                splitdict['D'].append(teminfo1+'\t'+ teminfo2+ '\t'+ teminfo3)
             




def main():
    args = sys.argv
    if len(args) < 3:
        usage(args)

    #pdb.set_trace()

    juncdict, gtfdict = {}, {}
    ReadGtf(args[1], gtfdict)
    ReadJunc(args[2], juncdict)
    juncdict = FiltJuncdic(gtfdict, juncdict)

    splitdict = {'A':[], 'B':[], 'C':[], 'D':[]}
    for chrom in Chrom:
        if chrom in juncdict.keys():
            s_exonlist = GetExonList(gtfdict, chrom)
            e_exonlist = GetExonList(gtfdict, chrom, 1)
            junclist, junclen = juncdict[chrom], len(juncdict[chrom])
            for i in xrange(junclen):
                if junclist[i][2] == 'fwd':
                    ConditionAC(gtfdict[chrom], chrom, junclist[i], e_exonlist, splitdict)
                else:
                    try:
                        ConditionBD(gtfdict[chrom], chrom, (junclist[i], junclist[i+1]), s_exonlist, splitdict)
                    except IndexError:
                        pass


    print 'Start to write the splitinfo ...'
    for k in splitdict:
        outfile = k + '.junc'
        fp = open(outfile, 'w')
        for e in splitdict[k]:
            fp.write( "%s\n" %(e) )
        fp.close()




if __name__ == '__main__':
    main()
