"""
  Function:  parse
  --------------------
  assembly of functions to parse text files
  
  momo.sander@googlemail.com
"""                              
def col2fltlist(path, col, header): 
  i = 0
  if header == True:
    i =1
  infile  = open(path, 'r')
  lines = infile.readlines()
  ll = []
  for line in lines[i:]:
    elements = line.split('\t')
    element = elements[col].rstrip('\n')
    ll.append(float(element))
  return ll  


def col2intlist(path ,sep, header, idx): 
  i = 0
  if header == True:
    i =1
  infile  = open(path, 'r')
  lines = infile.readlines()
  ll = []
  for line in lines[i:]:
    elements = line.split(sep)
    element = elements[idx].rstrip('\n')
    ll.append(int(element))
  return ll 




def col2keys(path,col, header): 
  "use cols2dict"


def col2list(path ,sep, header, idx): 
  i = 0
  if header == True:
    i =1
  infile  = open(path, 'r')
  lines = infile.readlines()
  ll = []
  for line in lines[i:]:
    elements = line.split(sep)
    element = elements[idx].rstrip()
    ll.append(element)
  return ll                                                                     


def rdstatLogs(path): 

  infile = open(path, 'r')
  lines = infile.readlines()
  elements = lines[0].split('\t') 
  al = float(elements[2])
  minx = int(elements[4])
  return(al, minx)



def cols2tup(path, header, (idxTup)):
  tups = []
  i = 0
  if header == True:
    i =1
  infile = open(path, 'r')
  lines = infile.readlines()

  for line in lines[i:]:
    elements = line.split('\t')
    tmp = []
    for idx in idxTup:
      key = elements[idx].rstrip('\n')
      tmp.append(key)
    tups.append(tuple(tmp))
  return tups



def cols2dict(path, symbol, header, keyIndex, valIndex):
  dctn = {}
  i = 0
  if header == True:
    i =1
  infile = open(path, 'r')
  lines = infile.readlines()

  for line in lines[i:]:
    elements = line.split(symbol)
    key = elements[keyIndex].strip('\n\"')
    value = elements[valIndex].strip('\n\"')
    try:
      value = int(elements[valIndex])
    except ValueError:
      pass
    dctn[key] = value
  return dctn

