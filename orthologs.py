"""
    Function: orthologs

    Constructs a distribution of differences in small molecule binding
    between pairs of orthologs (for now: human-to-rat). 
    -------------------

    Felix Kruger
    momo.sander@googlemail.com
"""

def orthologs():
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
    comparaOrthologs = "data/orthologs_%s.txt"% params['vsCompara']
    comparaHumanIds = "data/humanIds_%s.txt"% params['vsCompara']
    comparaRatIds = "data/ratIds_%s.txt"% params['vsCompara']
    # Assign output filenames.
    dictFile = "data/ortho_compDict_%s.pkl"% params['release']
    results = "data/orthologs_%s_%s.tab"%(params['release'], params['vsCompara'])
    orthoTab = "data/orthologTable_%s.txt" % params['vsCompara']
    # Create output files.
    humanLkp = readIds.readIds(comparaHumanIds)
    ratLkp = readIds.readIds(comparaRatIds)
    mkHomologTable.homologTable(comparaOrthologs, orthoTab, humanLkp, ratLkp) 
    query = queries.paralogs(orthoTab)
    acts= queryDevice.queryDevice(query, params)
    mkDict.activities(acts, dictFile)
    writePairs.homologMedian(params['homologyTypeOrthologs'], orthoTab, dictFile, results)
    # Annotate output files.
    align.pfam_a(results, params)
    align.bSite(results, params)
    addProperties.addMolweight('molregno', results, params)
    addProperties.addTargetClass("L1","accession1", results, params)
    addProperties.addPrefName("accession1", results, params)
    addProperties.addPrefName("accession2", results, params)

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 1: # the program name and the two arguments 
        sys.exit("All parameters are passed through gla.yaml")
    orthologs()



