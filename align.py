"""
    Function:  alignSeqs
    --------------------
For each pair in the paralogTable identifies the bindingSite containing domain 
from the binding_sites table of chembl and then ligns them using the EMBOSS
needle algorithm.
"""          
def fullSeq(path, params):
    import pickle
    import queryDevice
    import needle
    import random
    inFile = open( path, 'r')
    lines = inFile.readlines()
    inFile.close()
    out = open(path ,'w')
    out.write("%s\tseq_id\tseq_sim\n"%lines[0].rstrip('\n'))
    seqIdDict = {}
    for line in lines[1:]:                
        elements = line.split('\t')
        proteinAcc_1 = elements[0]
        proteinAcc_2 = elements[1]
        pairName = ('_').join([proteinAcc_1, proteinAcc_2])
        try:            
            (seqSim, seqId) = seqIdDict[pairName]
            out.write("%s\t%s\t%s\n"%(line.rstrip('\n'), seqId, seqSim))
            continue
        except KeyError:
            #print "aligning sequences of: %s\t%s"%(proteinAcc_1, proteinAcc_2)
            pass
        data = queryDevice.queryDevice("SELECT td.protein_sequence, td.protein_accession FROM target_dictionary td WHERE td.protein_accession IN ('%s')"% "','".join([proteinAcc_1, proteinAcc_2]), params)
        lkp = {}
        for entry in data:
            lkp[entry[1]] = entry[0]
        try:
            seq_1 = lkp[proteinAcc_1]
            seq_2 = lkp[proteinAcc_2]     
        except KeyError:
            seqIdDict[pairName] = (None, None)
            out.write("%s\t%s\t%s\n"%(line.rstrip('\n'), None, None))
            continue
        ################################################
        # Align the sequences using needle from EMBOSS.
        needleReport = needle.needle(params['needlepath'], seq_1, seq_2)
        # Parse the output of the alignment
        (seqSim, seqId) = needle.parseNeedle(needleReport)
        seqIdDict[pairName] = (seqSim, seqId)
        out.write("%s\t%s\t%s\n"%(line.rstrip('\n'), seqId, seqSim))
    out.close()
  



def pfam_a(path, params):
    import queryDevice
    import needle
    inFile = open( path, 'r')
    lines = inFile.readlines()
    inFile.close()
    out = open(path ,'w')
    out.write("%s\tdom_seq_id\tdom_seq_sim\tpfam_1\tpfam_2\n"%lines[0].rstrip('\n'))
    seqIdDict = {}
    for line in lines[1:]:
        elements = line.split('\t')
        proteinAcc_1 = elements[0]
        proteinAcc_2 = elements[1]
        pairName = ('_').join([proteinAcc_1, proteinAcc_2])
        try:            
            (seqSim, seqId, pfam_1, pfam_2) = seqIdDict[pairName]
            out.write("%s\t%s\t%s\t%s\t%s\n"%(line.rstrip('\n'), seqId, seqSim, pfam_1, pfam_2))
            continue
        except KeyError:
            #print "aligning bs_containing domains of: %s\t%s"%(proteinAcc_1, proteinAcc_2)
            pass
        data_1 = queryDevice.queryDevice("""SELECT DISTINCT mp.pfam_a,
                             pd.start, pd.end, td.protein_sequence, td.protein_accession 
                         FROM map_pfam mp 
                         JOIN pfam_domains pd 
                         ON mp.pfam_a = pd.pfam_a 
                         JOIN target_dictionary td 
                         ON td.protein_accession = mp.protein_accession
                         WHERE mp.protein_accession = '%s' 
			 AND pd.protein_accession =  '%s' """ % (proteinAcc_1, proteinAcc_1), params)
        data_2 = queryDevice.queryDevice("""SELECT DISTINCT mp.pfam_a,
                             pd.start, pd.end, td.protein_sequence, td.protein_accession 
                         FROM map_pfam mp 
                         JOIN pfam_domains pd 
                         ON mp.pfam_a = pd.pfam_a 
                         JOIN target_dictionary td 
                         ON td.protein_accession = mp.protein_accession
                         WHERE mp.protein_accession = '%s' 
                         AND pd.protein_accession = '%s' """ % (proteinAcc_2, proteinAcc_2),params)

        lkp = {}
        for entry in data_1 + data_2:
            (pfam, start, end, fullSeq, acc) = (entry[0], entry[1], entry[2], entry[3], entry[4])
            seq = fullSeq[start:end]
            lkp[acc] = (seq, pfam)
        try:
            seq_1 = lkp[proteinAcc_1][0]
            pfam_1 =lkp[proteinAcc_1][1]
            seq_2 = lkp[proteinAcc_2][0]
            pfam_2 = lkp[proteinAcc_2][1]
        except KeyError:
            out.write("%s\t%s\t%s\t%s\t%s\n"%(line.rstrip('\n'), None, None, None, None))
            continue
        ################################################
        # Align the sequences using needle from EMBOSS.
        needleReport = needle.needle(params['needlepath'], seq_1, seq_2)
        ################################################
        # Parse the output of the alignment
        (seqSim, seqId) = needle.parseNeedle(needleReport)
        seqIdDict[pairName] = (seqSim, seqId, pfam_1, pfam_2)
        out.write("%s\t%s\t%s\t%s\t%s\n"%(line.rstrip('\n'), seqId, seqSim, pfam_1, pfam_2))
    out.close()




def sliceBS(proteinAcc, proteinAli, sitePositions):
    bsDict = {}
    for acc in proteinAcc.keys():
        dom = proteinAcc[acc]
        seq = proteinAli[dom]
        bsSeq = ''
        for pos in sitePositions:
            bsSeq = bsSeq + seq[pos]
        bsDict[acc] = bsSeq
    return bsDict


def parseAcc(ksAccFile):
    dctn = {}
    infile = open(ksAccFile, 'r')
    lines = infile.readlines()
    for line in lines[1:]:
        elements = line.split('\t')
        accs = elements[11].strip()
        value = int(elements[0])
        for key in accs.split(', '):
            dctn[key] = value
    return dctn

def parseGSAli(aliFile):
    import re
    inFile = open(aliFile,'r')
    lines = inFile.readlines()
    out = open('data/reformatAli.log', 'w')
    alignDict = {}
    for line in lines:
        # Check if line holds identifier.
        if re.match('>',line):
            domId = int(line.split('_')[2].rstrip('\n'))
            seq = ''
            alignDict[domId]  = ''
        # Check if line holds sequence end.  
        elif re.search('\*',line):
            seq = ''.join([seq, line.rstrip('*\n')])
            alignDict[domId] = seq
        # Process any other line.                                
        else:
            seq = ''.join([seq, line.rstrip('\n')])
    for dom in alignDict.keys():
        out.write('%s\t%s\n'% (domId, alignDict[domId]))
    out.close()
    return alignDict 

def parseKSAli(aliFile):
    inFile = open(aliFile,'r')
    lines = inFile.readlines()
    alignDict = {}
    for line in lines[1:]:
        elements = line.split('\t')
        displayName = elements[0]
        domId = int(displayName.split('_')[-1])
        aln = elements[1].rstrip()
        alignDict[domId] = aln
    return alignDict 

def bSite(path, params):
    """Aligns binding site residues of GPCRs and kinases."""
    import queryDevice
    import parse
    import needle
    ksAccFile = params['ksAccFile']
    gsAccFile = params['gsAccFile'] 
    ksAliFile = params['ksAliFile'] 
    gsAliFile = params['gsAliFile']
    sitePositions = params['sitePositions']
    kinaseAcc = parseAcc(ksAccFile)
    gpcrAcc = parseAcc(gsAccFile)
    kinaseAli = parseKSAli(ksAliFile)
    gpcrAli = parseGSAli(gsAliFile)
    kinaseDict = sliceBS(kinaseAcc, kinaseAli, sitePositions['kinase'])
    gpcrDict = sliceBS(gpcrAcc, gpcrAli, sitePositions['gpcr'])
    # Processing the table of homologous pairs.
    inFile = open( path, 'r')
    lines = inFile.readlines()
    inFile.close()
    out = open(path ,'w')
    out.write("%s\tbs_seq_id\tbs_seq_sim\n"%lines[0].rstrip('\n'))
    seqIdDict = {}
    for line in lines[1:]:
        elements = line.split('\t')
        proteinAcc_1 = elements[0]
        proteinAcc_2 = elements[1]
        pairName = ('_').join([proteinAcc_1, proteinAcc_2])
        try:
            (seqSim, seqId) = seqIdDict[pairName]
            out.write("%s\t%s\t%s\n"%(line.rstrip('\n'), seqId, seqSim))
            continue
        except KeyError:
            pass    
        if proteinAcc_1 in kinaseDict or proteinAcc_2 in kinaseDict:
            (seq_1, seq_2) = (kinaseDict[proteinAcc_1], kinaseDict[proteinAcc_2])
        elif proteinAcc_1 in gpcrDict and proteinAcc_2 in gpcrDict:
            (seq_1, seq_2) = (gpcrDict[proteinAcc_1], gpcrDict[proteinAcc_2])
        else: 
            out.write("%s\t%s\t%s\n"%(line.rstrip('\n'), None, None))
            continue
        ################################################
        # Align the sequences using needle from EMBOSS.
        needleReport = needle.ungapped_needle(params['needlepath'], seq_1, seq_2)
        ################################################
        (seqSim, seqId) = needle.parseNeedle(needleReport)
        seqIdDict[pairName]=(seqSim, seqId)
        out.write("%s\t%s\t%s\n"%(line.rstrip('\n'), seqId, seqSim))
    out.close()



