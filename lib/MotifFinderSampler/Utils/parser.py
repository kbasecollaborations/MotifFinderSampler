import sys                                                                                                      
import os                                                                                                       
import json                            
import re
import numpy as np
from Bio import motifs
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from io import StringIO

import shutil
import subprocess

class SamplerUtil:
  def __init__(self):
      pass

  def parse_sampler_output(self, path):
      outputFileList = []
 
      seqflag=False
      pwmDict={}
      a=[]
      c=[]
      g=[]
      t=[]
      for filename in os.listdir(path):
          outputFileList.append(path + '/' + filename)
          if(filename=="SeqSet.matrix"):
             outputFilePath=path+'/'+filename
             print(outputFilePath)
             samplerFile = open(outputFilePath,'r')
             seqid=''
             for line in samplerFile:
                 line=line.rstrip()
                 if(line.startswith("#ID")):
                   out=line.split()
                   
                   seqid=out[2]
                   
                   if(len(pwmDict) !=0 ):
                      print(pwmDict)

                   seqflag=True
                  
                   seq=line.split("\t")
              
                   pwmDict['A']=a
                   pwmDict['C']=c
                   pwmDict['G']=g
                   pwmDict['T']=t
                   if(len(pwmDict)):
                      print(seqid)
                      print(pwmDict)
                   a=[]
                   c=[]
                   g=[]
                   t=[]
                   pwmDict={}
                 if(seqflag):
                     if((line.strip() != '') and (line[0]).isdigit()):
                        out=line.split()
                        a.append(float(out[0]))
                        c.append(float(out[1]))
                        g.append(float(out[2]))
                        t.append(float(out[3]))
             else:
                pwmDict['A']=a
                pwmDict['C']=c
                pwmDict['G']=g
                pwmDict['T']=t
                print(seqid)
                print(pwmDict)
                a=[]
                c=[]
                g=[]
                t=[]
                pwmDict={}    
      return pwmDict



Su=SamplerUtil()
out=Su.parse_sampler_output("/home/manish/Desktop/reorganization/MotifFinderSampler/test_local/workdir/tmp/sampler_out")
#print(out)

