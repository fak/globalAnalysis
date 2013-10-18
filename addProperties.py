"""
    Function:  addMolweight

    adds the column Molweight to paralogFull.tab and writes out to 
    paralogMolweight.tab
    --------------------

    Felix Kruger
    momo.sander@googlemail.com
"""



def addMolweight(key, path, params):
    import queryDevice
    import os
    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()
    out = open('_'.join([path,"sed"]) ,'w')
    out.write('%s\tmolweight\n'%lines[0].rstrip('\n') )
    header = lines[0].split('\t')
    molregnos = {}
    for i, col in enumerate(header):
        if col == key:
            idx = i
            break
    for line in lines[1:]:
        elements = line.split('\t')
        molregno = int(elements[idx])
        molregnos[molregno] = 0 
    molstr = "','".join(map(str, molregnos.keys()))
    print "Looking up mw_freebase for ", len(molregnos.keys()), "cmpds."
    data = queryDevice.queryDevice("SELECT distinct molregno, mw_freebase FROM compound_properties WHERE molregno IN('%s')"% molstr, params)
    for tup in data:
        molregno = int(tup[0])
        molweight = tup[1]
        molregnos[molregno] = molweight
    for line in lines[1:]:
        elements = line.split('\t')
        molregno = int(elements[idx])
        molweight = molregnos[molregno]
        out.write("%s\t%s\n"%(line.rstrip('\n'), molweight ))
    out.close()
    os.system('mv %s %s'% ('_'.join([path,"sed"]), path))


def addChembl_id(key, path, params):
    import queryDevice
    import os
    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()
    out = open('_'.join([path,"sed"]) ,'w')
    out.write('%s\tchembl_id\n'%lines[0].rstrip('\n') )
    header = lines[0].split('\t')
    molregnos = {}
    for i, col in enumerate(header):
        if col == key:
            idx = i
            break
    for line in lines[1:]:
        elements = line.split('\t')
        molregno = int(elements[idx])
        molregnos[molregno] = 0
    molstr = "','".join(map(str, molregnos.keys()))
    print "Looking up chembl_id for ", len(molregnos.keys()), "cmpds."
    data = queryDevice.queryDevice("SELECT distinct molregno, chembl_id FROM molecule_dictionary WHERE molregno IN('%s')"% molstr, params)
    for tup in data:
        molregno = int(tup[0])
        chembl_id = tup[1]
        molregnos[molregno] = chembl_id
    for line in lines[1:]:
        elements = line.split('\t')
        molregno = int(elements[idx])
        chembl_id = molregnos[molregno]
        out.write("%s\t%s\n"%(line.rstrip('\n'), chembl_id ))
    out.close()
    os.system('mv %s %s'% ('_'.join([path,"sed"]), path))



def addTargetClass(level, key, path, params):
    import queryDevice
    import os
    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()
    out = open('_'.join([path,"sed"]) ,'w')
    out.write('%s\ttarget_class_%s\n'%(lines[0].rstrip('\n'), level))
    header = lines[0].split('\t')
    accessions = {}
    for i, col in enumerate(header):
        if col == key:
            idx = i
            break
    for line in lines[1:]:
        elements = line.split('\t')
        acc = elements[idx]
        accessions[acc] = 0
    accStr = "','".join(map(str, accessions.keys()))
    data = queryDevice.queryDevice("SELECT protein_accession, %s FROM target_class tc JOIN target_dictionary td ON td.tid = tc.tid  WHERE protein_accession IN('%s') "%(level, accStr), params)
    for tup in data:
        acc = tup[0]
        targetClass = tup[1]
        accessions[acc] = targetClass
    for line in lines[1:]:
        elements = line.split('\t')
        acc = elements[idx]
        targetClass = accessions[acc]
        out.write("%s\t%s\n"%(line.rstrip('\n'), targetClass ))
    out.close()
    os.system('mv %s %s'% ('_'.join([path,"sed"]), path))


def addSeq100(path):
    import os
    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()
    out = open('_'.join([path,"sed"]) ,'w')
    out.write('%s\tseqId\n'%lines[0].rstrip('\n'))
    for line in lines[1:]:   
        out.write("%s\t%s\n"%(line.rstrip('\n'), 100 ))    
    out.close()
    os.system('mv  %s %s'% ('_'.join([path,"sed"]), path))



def addPrefName(key, path, params):
    import queryDevice
    import os
    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()
    out = open('_'.join([path,"sed"]) ,'w')
    out.write('%s\tprefName_%s\n'%(lines[0].rstrip('\n'), key))
    header = lines[0].split('\t')
    accessions = {}
    for i, col in enumerate(header):
        if col == key:
            idx = i
            break
    for line in lines[1:]:
        elements = line.split('\t')
        acc = elements[idx]
        accessions[acc] = 0
    accStr = "','".join(map(str, accessions.keys()))
    data = queryDevice.queryDevice("SELECT protein_accession, pref_name FROM target_dictionary WHERE protein_accession IN('%s') "% accStr, params)
    for tup in data:
        acc = tup[0]
        prefName = tup[1]
        accessions[acc] = prefName

    for line in lines[1:]:
        elements = line.split('\t')
        acc = elements[idx]
        prefName = accessions[acc]
        out.write("%s\t\"%s\"\n"%(line.rstrip('\n'), prefName))
    out.close()
    os.system('mv %s %s'% ('_'.join([path,"sed"]), path))



