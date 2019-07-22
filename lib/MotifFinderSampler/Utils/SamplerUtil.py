import sys                                                                                                      
import os                                                                                                       
import json                            
import re
import numpy as np
from Bio import motifs
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from io import StringIO
from installed_clients.DataFileUtilClient import DataFileUtil
import shutil
import subprocess

class SamplerUtil:
  def __init__(self):
      pass

  def build_sampler_motif_command(self, inputFilePath, motiflen, prb):
      cmd1 = 'cp -r /kb/module/deps/kb_sampler/CreateBackgroundModel /kb/module/work/tmp/sampler_out'
      cmd2 = 'cp -r /kb/module/deps/kb_sampler/MotifSampler /kb/module/work/tmp/sampler_out'
      subprocess.call('mkdir /kb/module/work/tmp/sampler_out', shell=True)
      subprocess.call(cmd1, shell=True)
      subprocess.call(cmd2, shell=True)
      os.chdir("/kb/module/work/tmp/sampler_out")
      os.system("/kb/module/work/tmp/sampler_out/CreateBackgroundModel -f /kb/module/work/tmp/SeqSet.fa -b /kb/module/work/tmp/sampler_out/SeqSet.bg")
      command = '/kb/module/work/tmp/sampler_out/MotifSampler -f /kb/module/work/tmp/SeqSet.fa -b /kb/module/work/tmp/sampler_out/SeqSet.bg -o /kb/module/work/tmp/sampler_out/SeqSet.out -m /kb/module/work/tmp/sampler_out/SeqSet.matrix -w 8'
      return command

  def run_sampler_command(self, command):
      os.system(command)

  def write_obj_ref(self, path, obj_ref):
      file = open(path+"/sampler_obj.txt","w")
      file.write(obj_ref)
      file.close() 

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

      
      motifDict['PWM'] = []
      motifDict['PFM'] = []

      
      pwmList=[]
      pwmDict={}
      pfmDict={}
      rowList = []
      rowDict={}
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

                   motifDict['Motif_Locations'] = []
                   locList=[]
                   seqflag=True
                   seq=line.split("\t")
                   consensus=(seq[1]).split(" ")[1]
                   pwmDict['A']=[]
                   pwmDict['C']=[]
                   pwmDict['G']=[]
                   pwmDict['T']=[]

                   pfmDict['A']=[]
                   pfmDict['C']=[]
                   pfmDict['G']=[]
                   pfmDict['T']=[]
                   motifDict['PWM']=pwmDict
                   motifDict['PFM']=pfmDict
                   motifDict['Iupac_sequence']=consensus
                   #print(consensus)
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

  def UploadFromSampler(self, callback_url, params):
          """
          :param params: instance of type "UploadmfmdInParams" -> structure:
             parameter "path" of String, parameter "ws_name" of String,
             parameter "obj_name" of String
          :returns: instance of type "UploadOutput" -> structure: parameter
             "obj_ref" of String
          """
          # ctx is the context object
          # return variables are: output
          #BEGIN UploadFromSampler
          print('Extracting motifs')
          motifList = self.parse_sampler_output(params['path'])
          print(motifList)
       
          MSO = {}
          MSO=motifList
        
          dfu = DataFileUtil(callback_url)
          save_objects_params = {}
          save_objects_params['id'] = dfu.ws_name_to_id(params['ws_name'])
          save_objects_params['objects'] = [{'type': 'KBaseGeneRegulation.MotifSet' , 'data' : MSO , 'name' : params['obj_name']}]

          info = dfu.save_objects(save_objects_params)[0]
          print('SAVED OBJECT')
          print(info)
          motif_set_ref = "%s/%s/%s" % (info[6], info[0], info[4])
          print(motif_set_ref)
          output = {'obj_ref' : motif_set_ref}
          print(output)

        
          #exit("test")
          #END UploadFromSampler

          # At some point might do deeper type checking...
          if not isinstance(output, dict):
              raise ValueError('Method UploadFrommfmd return value ' +
                             'output is not type dict as required.')
          # return the results
          return [output]

#Su=SamplerUtil()
#out=Su.parse_sampler_output("/home/manish/Desktop/reorganization/MotifFinderSampler/test_local/workdir/tmp/sampler_out")
#print(out)

