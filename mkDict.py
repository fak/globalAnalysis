"""
    Function:  mkDict
    Converts input from the SQL statement into a dictionary that allows for the
    comparison of potencies measured in different species. Where necessary, 
    potencies are converted to negative-logarithms. 
    --------------------

    Felix Kruger
    momo.sander@googlemail.com
"""

def activities(acts, dictFile):  
    import pickle                                                       
    import math
    import numpy as np 
    import random
    compDict = {}
    lookUp= {}
    molregs = list(acts) 
    print "found %s activities, processing" % len(molregs) 
    random.shuffle(molregs)
    for molreg in molregs:
        pAfnty=0  
        molregno = int(molreg[0])
        standardValue = molreg[1]
        try:
            standardValue = float(standardValue)
        except:
            continue
        standardValue = float(standardValue) 
        if standardValue <=0:  
            continue
        standardType = molreg[2]
        standardUnit = molreg[3]
        assayType = molreg[4]
        assayId = molreg[5]                   
        prefName = molreg[6]
        uniprot = molreg[7]
        reference = str(int(molreg[8]))
        if standardType in ['IC50', 'EC50', 'Ki'] and standardUnit == 'nM':
            pAfnty = (-1) * math.log10(standardValue/float(1000000000))
        elif standardType in ['pIC50', 'pEC50','pKi']:
            pAfnty = standardValue         
        else:
            continue
        if standardType in ['Ki', 'pKi']:
            pAfnty = pAfnty - 0.355
        if molregno not in compDict.keys():
            compDict[molregno] = {}
        if uniprot not in compDict[molregno].keys():
            compDict[molregno][uniprot]={}     
            compDict[molregno][uniprot]['pAfnty'] = []
            compDict[molregno][uniprot]['references']=[]
            compDict[molregno][uniprot]['assayId']=[]
            compDict[molregno][uniprot]['molregno']=[]
            compDict[molregno][uniprot]['prefName'] = []      
        compDict[molregno][uniprot]['pAfnty'].append(pAfnty)
        compDict[molregno][uniprot]['references'].append(reference)
        compDict[molregno][uniprot]['assayId'].append(assayId)
        compDict[molregno][uniprot]['molregno'].append(molregno)
        compDict[molregno][uniprot]['prefName'].append(prefName) 
    out = open(dictFile, 'w')
    pickle.dump(compDict, out)  
    out.close()



