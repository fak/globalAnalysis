"""
    Function: needle 
    for each line from the inputFile, determines the sequence identity of the 
    biniding site containing domains in the pair of orthologs. Results are written
    to a table. (path = originalFile + 'Seq.tab')

    --------------------
""" 

def needle(needlepath, seq1, seq2):
    import subprocess
    process = subprocess.Popen(['%s'%needlepath , '-gapopen=10','-datafile=EBLOSUM62','asis:%s'%seq1,'asis:%s'%seq2, \
                '-gapextend=0.5','outfile=stdout'], shell=False, stdout=subprocess.PIPE, stdin=subprocess.PIPE)  
    return process.communicate()[0]
if __name__ == '__main__':
    needle()


def ungapped_needle(needlepath, seq1, seq2):
    import subprocess
    process = subprocess.Popen(['%s'%needlepath , '-gapopen=100','-datafile=EBLOSUM62','asis:%s'%seq1,'asis:%s'%seq2, \
                '-gapextend=10','outfile=stdout'], shell=False, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    return process.communicate()[0]
if __name__ == '__main__':
    needle()


def parseNeedle(needleReport):
    """Extract sequence identity and sequence similarity from needle report."""
    import re
    seqSim = re.search('Similarity.*', needleReport)
    seqSim = seqSim.group(0)
    seqSim = seqSim.split('(')[1]
    seqSim = float(seqSim.split('%')[0])
    seqId = re.search('Identity.*', needleReport)
    seqId = seqId.group(0)
    seqId = seqId.split('(')[1]
    seqId = float(seqId.split('%')[0])
    return seqSim, seqId

