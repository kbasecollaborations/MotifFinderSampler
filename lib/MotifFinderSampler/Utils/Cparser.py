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

  def parse_sampler_output(self, path):
      outputFileList = []
 
      seqflag=False
      motifList={}
      motifDict={}
      locList=[]
      alphabet=['A','C','G','T']
      print(path)
  
      motifSet=[]
      motifList['Condition']='temp'
      motifList['SequenceSet_ref']='123'

      background={}
      background['A']=0.0
      background['C']=0.0
      background['G']=0.0
      background['T']=0.0

      
      
      
      pwmList=[]
      
      pfmDict={}
      rowList = []
      rowDict={}
      
      matrix=self.parse_matrix_output(path)

      for filename in os.listdir(path):
          outputFileList.append(path + '/' + filename)
          if(filename=="SeqSet.out"):
             outputFilePath=path+'/'+filename
             print(outputFilePath)
             samplerFile = open(outputFilePath,'r')
             for line in samplerFile:
                 line=line.rstrip()
                 if(line.startswith("#id:")):
                   #print(motifDict)
                   if(len(motifDict) !=0 ):
                      motifSet.append(motifDict)
                      pwmDict={}
                      motifDict={}
                      

                   motifDict['Motif_Locations'] = []
                   locList=[]
                   seqflag=True
                   seq=line.split("\t")
                   consensus=(seq[1]).split(" ")[1]
                   seqid=(seq[0]).split(" ")[1]
                  
                   pwmDict={}
                   pwmDict['A']=[]
                   pwmDict['C']=[]
                   pwmDict['G']=[]
                   pwmDict['T']=[]
                   motifDict['PWM']=self.get_pwm(matrix[seqid])
                   print(seqid)
                   print(self.get_pwm(matrix[seqid]))
                   motifDict['PFM']=pfmDict
                   motifDict['Iupac_sequence']=consensus
                   
                 if(seqflag):
                     if(line.endswith(";")):
                        out=line.split(" ")
                        rec=(out[0]).split("\t")
                        seqid=rec[0]
                        orientation=rec[6]
                        if(orientation == "+"):
                           seq_start=rec[3]
                           seq_end=rec[4]
                        else:
                           seq_start=rec[4]
                           seq_end=rec[3]
                        sequence=(out[3]).replace(";", "").replace("\"", "")
                        locDict={}
                        locDict['sequence_id']=seqid;
                        locDict['start']=int(seq_start);
                        locDict['end']=int(seq_end);
                        locDict['sequence']=sequence;
                        locDict['orientation']=orientation;
                        motifDict['Motif_Locations'].append(locDict)
                        #print(seqid+"\t"+sequence+"\t"+strand+"\t"+start+"\t"+stop)
                   
                        locList.append(line)
             else:
                 #print(motifDict)
                 motifSet.append(motifDict)

      motifList['Motifs']=motifSet
      motifList['Background']=background
      motifList['Alphabet']=alphabet  
      #print(motifList)     
      return motifList

                
Su=SamplerUtil()
#out=Su.parse_matrix_output("/home/manish/Desktop/reorganization/MotifFinderSampler/test_local/workdir/tmp/sampler_out")
#dict=Su.get_pwm(out['box_1_1_CCTTCnTC'])
out=Su.parse_sampler_output("/home/manish/Desktop/reorganization/MotifFinderSampler/test_local/workdir/tmp/sampler_out")
print(out)

