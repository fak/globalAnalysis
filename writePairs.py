"""
    Module:  writePairs

    Writes out pairs of potency measurements for comparison i) between assays
    ii) between orthologs, iii) between paralogs.
    --------------------

    Felix Kruger
    momo.sander@googlemail.com
"""


def homologMedian(homologyTypes, homologTable, dictFile, outfile):
    import pickle
    import numpy as np
    infile = open(homologTable, 'r')
    lines = infile.readlines()
    infile.close()
    infile = open(dictFile , 'r')
    compDict = pickle.load(infile)
    infile.close()
    out = open(outfile ,'w')
    out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%("accession1","accession2", "seq_id",\
                                               "molregno","afnty1", "afnty2","ref1","ref2"))
    for line in lines[1:]:
        elements = line.split('\t')
        seqId = elements[4]
        viceversaSeqId = elements[5].rstrip('\n')
        if viceversaSeqId > seqId:
            seqId = viceversaSeqId
        try:
            seqId = float(seqId)
        except:
            print "Couldn't get seq Id: ", seqId
            continue
        if not elements[2] in homologyTypes:
            print elements[2]
            continue
        if seqId<=0:
            continue
        accession1 = elements[0]
        accession2 = elements[1]
        for molregno in compDict.keys():
            if accession1 in compDict[molregno].keys() and  \
            accession2 in compDict[molregno].keys():
                ref1 = compDict[molregno][accession1]['references']
                ref1 = ','.join(ref1)
                ref2 = compDict[molregno][accession2]['references']
                ref2 = ','.join(ref2)
                pAfnty1 = np.median(compDict[molregno][accession1]['pAfnty'])
                pAfnty2 = np.median(compDict[molregno][accession2]['pAfnty'])
                #print "writing data for pair: %s\t%s"%(accession1, accession2)
                out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"%( accession1,accession2, \
                                                        seqId, molregno, pAfnty1, pAfnty2, ref1, ref2))
    out.close()



def interAssaySampled(path, dictFile):
    import pickle
    import random
    pick = open(dictFile , 'r')
    compDict = pickle.load(pick)
    pick.close()
    out = open( path,'w')
    out.write("%s\t%s\t%s\t%s\t%s\t%s\n"%("prefName","accession","afnty1",\
       "afnty2", "molregno", "measured"))
    for molregno in compDict.keys():
        for uniprot in compDict[molregno].keys():
            if len(compDict[molregno][uniprot]['pAfnty'])>1:
                groupSize = len(compDict[molregno][uniprot]['pAfnty'])
                afnty = []
                for i in range(0,groupSize):
                    failCount = 0
                    for j in range(i+1, groupSize):
                        if compDict[molregno][uniprot]['assayId'][i] == compDict[molregno][uniprot]['assayId'][j]\
                        or compDict[molregno][uniprot]['pAfnty'][i] == compDict[molregno][uniprot]['pAfnty'][j]:
                            failCount += 1
                    if failCount ==0:
                        afnty.append(compDict[molregno][uniprot]['pAfnty'][i])
                groupSize = len(afnty)
                if groupSize > 1:
                    (afnty1,afnty2) = random.sample(afnty,2)
                    prefName=compDict[molregno][uniprot]['prefName'][1]
                    measured = len(afnty)
                    out.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(prefName, uniprot, afnty1,\
                                                    afnty2,  molregno,measured))
    out.close()

#### Function below was used in previous approach but not recommended.

def interAssayMedian(path, dictFile):
    import pickle
    import numpy as np
    import random
    pick = open(dictFile , 'r')
    compDict = pickle.load(pick)
    pick.close()
    out = open(path, 'w')
    out.write("%s\t%s\t%s\t%s\t%s\t%s\n"%("prefName","accession","afnty1", "afnty2", \
          "molregno", "measured"))
    for molregno in compDict.keys():
        for uniprot in compDict[molregno].keys():
            if len(compDict[molregno][uniprot]['pAfnty'])>1:
                groupSize = len(compDict[molregno][uniprot]['pAfnty'])
                afnty = []
                for i in range(0,groupSize):
                    failCount = 0
                    for j in range(i+1, groupSize):
                        if compDict[molregno][uniprot]['assayId'][i] == compDict[molregno][uniprot]['assayId'][j]\
                        or compDict[molregno][uniprot]['pAfnty'][i] == compDict[molregno][uniprot]['pAfnty'][j]:
                            failCount += 1
                    if failCount ==0:  # append pAfnty[i] only if there is no pAfnty[j] == pAfnty[i]
                        afnty.append(compDict[molregno][uniprot]['pAfnty'][i])
                groupSize = len(afnty)
                if groupSize > 1:
                    random.shuffle(afnty)
                    halfSize = groupSize/2 #deliberate rounding down
                    afnty1 = np.median(afnty[0:halfSize])
                    afnty2 = np.median(afnty[halfSize:groupSize])
                    prefName=compDict[molregno][uniprot]['prefName'][1]
                    measured = len(afnty)
                    out.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(prefName, uniprot, afnty1,\
                                                    afnty2,  molregno,measured))
    out.close()
