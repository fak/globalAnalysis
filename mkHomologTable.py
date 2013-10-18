
"""
    Function:  createParalogTable
    
    This function processes the output of an Ensembl query. It filters out
    paralogous pairs that have no Uniprot identifiers. It also removes pairs 
    that are duplicated through inversion. 
    --------------------
    Felix Kruger
    momo.sander@googlemail.com
"""

def createParalogTable(infile, outfile, humanIds, paralogIds):
    import random
    infile = open(infile, 'r')
    lines = infile.readlines()
    infile.close()
    out=open(outfile, 'w')
    out.write("humanId_1\thumanId_2\tancestor\ttype\tseqId\tviceversaSeqId\n")
    sims={}
    inverseList=[]
    lines = lines[1:] 
    random.shuffle(lines) # Do this to avoid alphabetic/sorting bias in the paralog random gouping.
    for line in lines:
        elements = line.split('\t')
        geneId1 = elements[0]
        geneId2 = elements[1]
        ancestor = elements[3]
        type = elements[2]
        seqId = elements[4]
        vvseqId = elements[5]
        if geneId1 in humanIds.keys() and geneId2 in paralogIds.keys(): 
            string = '_'.join([humanIds[geneId1],paralogIds[geneId2]]) 
            inverseString = '_'.join([humanIds[geneId2],paralogIds[geneId1]]) 
            if inverseString in inverseList:
                continue
            out.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(humanIds[geneId1],     \
               paralogIds[geneId2], ancestor, type, seqId, \
               vvseqId.rstrip('\n')))
            inverseList.append(string)                        
    out.close()



def homologTable(infile, outfile, ids, homoIds):
    import random
    infile = open(infile, 'r')
    lines = infile.readlines()
    infile.close()
    out=open(outfile, 'w')
    out.write("accession\thomolog_accession\tancestor\ttype\tseqId\tviceversaSeqId\n")
    lkp = {}
    lines = lines[1:]
    random.shuffle(lines) # Do this to avoid alphabetic/sorting bias in the paralog random grouping.
    for line in lines:
        elements = line.split('\t')
        geneId_1 = elements[0]  
        geneId_2 = elements[1]
        type = elements[2]
        ancestor = elements[3]
        seqId = elements[4]
        vvseqId = elements[5]
        if geneId_1 in ids.keys() and geneId_2 in homoIds.keys():  
            accessionStr = '_'.join(sorted([ids[geneId_1], homoIds[geneId_2]]))
            if accessionStr not in lkp.keys():
                lkp[accessionStr] = 0
                out.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(ids[geneId_1], \
                         homoIds[geneId_2], type, ancestor, seqId, vvseqId.rstrip('\n')))
    out.close()


