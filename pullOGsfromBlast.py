#!/usr/bin/env python
# usage: ./test.py genomepath blastoutputpath colrefOG colquery
# example: ./test.py /path/to/Hydra.fas /path/to/Hydra.blastout colrefOG colquery

#
#		example of tab-delimited blastoutput-- /path/to/Hydra.blastout:
#			NematostellaZZZ1  HydraZZZ3  123  456 
#			NematostellaZZZ2  HydraZZZ4  789  123
#
#			NematostellaZZZ1 is the reference OG (col 1); becomes name of fasta file
#			HydraZZZ3 (col 2) gets pulled from its genome and appended to NematostellaZZZ1.fas


import sys, os, subprocess
from os.path import basename
import re

def pullseqs(fastaseqheader, fastapath):
  p = subprocess.Popen(['selectSeqs.pl','-m',fastaseqheader,fastapath], stdout=subprocess.PIPE)
	# selectSeqs.pl -m Nematostella  $i >> Nemato_OGs.fas
  out, err = p.communicate()
  return(out)

def main():
  genomepath, blastoutput, colrefOG, colquery = sys.argv[1:]
  with open(blastoutput) as f1:  
    for i in f1.readlines():
      outfilename = i.split('\t')[int(colrefOG)-1] + '.fa'
      if os.path.exists(outfilename):
        apwr =  'a'
      else:
      	apwr = 'w'
      preseqid = i.split('\t')[int(colquery)-1]
      seqid = re.sub("\|","\\\|",preseqid) + '$' 
      myseq = pullseqs(seqid, genomepath)
      with open(outfilename,apwr) as myOGseqs:
        myOGseqs.write(myseq)

if __name__ == '__main__':
  main()
