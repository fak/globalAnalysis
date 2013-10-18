"""
  Function: readIds
  Reads the bioMart ID tables - it is important to use this rather than the 
  parse module because it will ignore empty cells in the uniprot column.
  --------------------

  Felix Kruger
  momo.sander@googlemail.com
"""

def readIds(inputFile):
  inFile = open(str(inputFile), 'r')
  lines = inFile.readlines()
  ids={}
  for line in lines[1:]:
    elements = line.split('\t')
    uniprot = elements[2].strip("\n")
    if len(uniprot) > 1:
      ids[elements[0]]=uniprot
  return ids

