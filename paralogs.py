"""
Function: paralogs
Constructs a distribution of differences in small molecule binding between pairs of human paralogs.
--------------------

Felix Kruger
momo.sander@googlemail.com
"""

def paralogs():
    import queryDevice
    import mkDict
    import queries
    import writePairs
    import os
    import align
    import addProperties
    import mkHomologTable
    import readIds
    import yaml
    # Read config file.
    paramFile = open('gla.yaml')
    params = yaml.safe_load(paramFile)
    needlepath = params['needlepath']
    vsCompara = params['vsCompara']
    release = params['release']
    comparaParalogs = "data/paralogs_%s.txt"% params['vsCompara']
    comparaHumanIds = "data/humanIds_%s.txt"% params['vsCompara']
    # Assign output filenames.
    dictFile = "data/para_compDict_%s.pkl"% params['release'] 
    results = "data/paralogs_%s_%s.tab"%(params['release'], params['vsCompara'])
    paraTab = "data/paralogTable_%s.txt" % vsCompara
    # Create output files.
    idLkp = readIds.readIds("data/humanIds_%s.txt"% vsCompara)
    mkHomologTable.homologTable(comparaParalogs, paraTab, idLkp, idLkp)
    query = queries.paralogs(paraTab)
    acts = queryDevice.queryDevice(query,params)
    mkDict.activities(acts, dictFile)
    writePairs.homologMedian(params['homologyTypeParalogs'], paraTab, dictFile, results)
    # Annotate output files.
    align.pfam_a(results, params)
    align.bSite(results, params)
    addProperties.addMolweight('molregno', results, params)
    addProperties.addTargetClass("L1","accession1", results, params)
    addProperties.addPrefName("accession1", results, params)
    addProperties.addPrefName("accession2", results, params)

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 1: # the program name and the two arguments i
        sys.exit("All parameters passed thorugh gla.yaml")
    paralogs()



