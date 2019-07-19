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
      #cwd=config['scratch']
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

  def parse_sampler_output(self, path):
      pfmList = []
      pfmDict={}
      outputFileList = []
      pfmMatrix=False
      seqflag=False
      motifList={}
      motifDict={}
      locList=[]
      alphabet=['A','C','G','T']
  
      '''motifSet=[]
      motifList['Condition']='temp'
      motifList['SequenceSet_ref']='123'

      background={}
      background['A']=0.0
      background['C']=0.0
      background['G']=0.0
      background['T']=0.0

      motifDict['Motif_Locations'] = []
      motifDict['PWM'] = []
      motifDict['PFM'] = []
      motiflen=0
      a=[]
      c=[]
      g=[]
      t=[]
      pwmList=[]
      pwmDict={}
      rowList = []
      rowDict={}

      for filename in os.listdir(path):
          outputFileList.append(path + '/' + filename)
          if(filename=="mfmd_out.txt"):
             outputFilePath=path+'/'+filename
             mfmdFile = open(outputFilePath,'r')
             for line in mfmdFile:
                 if(re.search("PPM Matrix",line)):
                    pfmMatrix=True
                 if(pfmMatrix):
                    if(line[0].isdigit()):
                       line=line.strip()
                       out=line.split()
                       pfmList.append(out)
                       a.append(out[0])
                       c.append(out[1])
                       g.append(out[2])
                       t.append(out[3])

                       rowList = []
                       rowList.append(('A',float(out[0])))
                       rowList.append(('C',float(out[1])))
                       rowList.append(('G',float(out[2])))
                       rowList.append(('T',float(out[3])))
                       rowDict['A']=float(out[0])
                       rowDict['C']=float(out[1])
                       rowDict['G']=float(out[2])
                       rowDict['T']=float(out[3])

                 if(re.search("PSSM Matrix",line)):
                    pfmMatrix=False
                 if(re.search("Sequences",line)):
                    seqflag=True
                 if(seqflag==True):
                    line=line.strip()
                    if(re.search('\*',line)):
                       seqflag=False
                    if((line) and not (line.startswith("Seq")) and not (line.startswith("*"))):
                       line=line.rstrip()
                       seq=line.split()
                       seqid=seq[0]
                       seq_start=int(seq[1])
                       seq_end=int(seq_start)+int(motiflen)
                       sequence=seq[2]
                       orientation='+'

                       locDict={}
                       locDict['sequence_id']=seqid;
                       locDict['start']=seq_start;
                       locDict['end']=seq_end;
                       locDict['sequence']=sequence;
                       locDict['orientation']=orientation;
                       motifDict['Motif_Locations'].append(locDict)

                  
                 if(re.search("Width",line)):
                    arr=line.split(" ")
                    motiflen=arr[1].split("\t")[0]

      a=[float(x) for x in a]
      c=[float(x) for x in c]
      g=[float(x) for x in g]
      t=[float(x) for x in t]
      pwmDict['A']=a
      pwmDict['C']=c
      pwmDict['G']=g
      pwmDict['T']=t

      pfmDict['A']=[]
      pfmDict['C']=[]
      pfmDict['G']=[]
      pfmDict['T']=[]

      motifStr = '>test\n'
      motifStr += 'A ' + str(a).replace(',','') + '\n'
      motifStr += 'C ' + str(c).replace(',','') + '\n'
      motifStr += 'G ' + str(g).replace(',','') + '\n'
      motifStr += 'T ' + str(t).replace(',','') + '\n'

      handle = StringIO(motifStr)

      BioMotif = motifs.read(handle, 'jaspar')
      motifDict['PWM']=pwmDict
      motifDict['PFM']=pfmDict
      motifDict['Iupac_sequence']=str(BioMotif.degenerate_consensus)
      motifSet.append(motifDict)                                            #keep in loop for multiple motifs

      motifList['Motifs']=motifSet
      motifList['Background']=background
      motifList['Alphabet']=alphabet'''

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
          #BEGIN UploadFrommfmd
          print('Extracting motifs')
          #motifList = MFU.parse_mfmd_output(params['path'])
          '''motifList = self.parse_mfmd_output(params['path'])
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
          print(output)'''

        
          #exit("test")
          #END UploadFrommfmd

          # At some point might do deeper type checking...
          if not isinstance(output, dict):
              raise ValueError('Method UploadFrommfmd return value ' +
                             'output is not type dict as required.')
          # return the results
          return [output]

