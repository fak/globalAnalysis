"""
    Function: interAssay

    This script generates a distribution of differencs observed between two experiments assessing the potency of identical compounds and targets
    --------------------

    Felix Kruger
    momo.sander@googlemail.com
"""
def interAssay():
    import queryDevice
    import mkDict
    import queries
    import writePairs
    import os
    import addProperties
    import readIds
    import string
    import yaml
    # Read config file.
    paramFile = open('gla.yaml')
    params = yaml.safe_load(paramFile)
    species = params['species']
    # Get information for all relevant activities from ChEMBL.
    for spec in species:
        specName =  string.replace(spec, ' ','_')
        dictFile =  "data/inter_compDict_%s_%s.pkl" % (specName, params['release'])
        results = "data/interAssay_%s_%s.tab" % (specName, params['release'])
        query = queries.activities(spec)
        acts = queryDevice.queryDevice(query, params) 
        mkDict.activities(acts, dictFile)
        writePairs.interAssaySampled(results, dictFile)
        addProperties.addMolweight("molregno", results, params)
        addProperties.addTargetClass("L1","accession", results, params)
        addProperties.addSeq100(results)


if __name__ == '__main__':

    import sys
    if len(sys.argv) != 1: # the program name and the two arguments i
        sys.exit("All parameters passed through gla.yaml")

    interAssay()

