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

  

  def is_eof(self, f):
    cur = f.tell()    # save current position
    f.seek(0, os.SEEK_END)
    end = f.tell()    # find the size of file
    f.seek(cur, os.SEEK_SET)
    return cur == end

  def get_pwm(self, lst):
      
      a=[]
      c=[] 
      g=[]
      t=[]
      pwmDict={}
      for line in lst:   
          out=line.split("\t")                  
          a.append(float(out[0]))
          c.append(float(out[1]))
          g.append(float(out[2]))
          t.append(float(out[3]))
                        
      pwmDict['A']=a
      pwmDict['C']=c
      pwmDict['G']=g
      pwmDict['T']=t
      
      return pwmDict

  def parse_matrix_output(self, path):
      data={}
      queryFile = open(path+'/'+"SeqSet.matrix",'r')
      qHeader=''
      qSequence=''
      while not self.is_eof(queryFile):
            qline = queryFile.readline()
             
            if (qline.startswith("#ID")):
                
                qHeader=qline.split(" ")[2]
                 
            else: 
                qSequence=qline  
                if((qline[0]).isdigit()):
                    qHeader=qHeader.strip()
                    qSequence=qSequence.strip() 
                    if(qHeader in data):
                       val=[]
                       val=data[qHeader]
                       val.append(qSequence)
                    else: 
                       val=[]
                       val.append(qSequence)
                       data[qHeader]=val

      return data
                
Su=SamplerUtil()
out=Su.parse_matrix_output("/home/manish/Desktop/reorganization/MotifFinderSampler/test_local/workdir/tmp/sampler_out")
dict=Su.get_pwm(out['box_1_1_CCTTCnTC'])
print(dict)

