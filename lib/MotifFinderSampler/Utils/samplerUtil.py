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

def build_mfmd_command(inputFilePath, motiflen, prb,config):
    
    shutil.copytree('/kb/module/deps/kb_sampler/', '/kb/module/work/tmp/sampler')
    print(prb)
    print(motiflen)
    cwd=config['scratch']
    os.chdir("/kb/module/work/tmp")
    command = 'java -jar mfmd/mfmd.jar ' + inputFilePath + ' ' + str(motiflen) + ' ' + str(prb)
    #command = 'java -jar /kb/module/work/tmp/mfmd/mfmd.jar ' + inputFilePath + ' ' + str(motiflen) + ' ' + str(prb)
    #print(command)
    return command

def run_mfmd_command(command):
    print(command)
    os.system('R CMD INSTALL /kb/module/deps/kb_mfmd/Rserve_1.7-3.1.tar.gz')
    os.system('/usr/lib/R/bin/Rscript /kb/module/deps/kb_mfmd/mfmd/script.R')
    os.system(command)

def parse_mfmd_output(path):
    pfmList = []
    pfmDict={}
    outputFileList = []
    pfmMatrix=False
    seqflag=False
    motifList={}
    motifDict={}
    locList=[]
    alphabet=['A','C','G','T']

    motifSet=[]
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
    motifList['Alphabet']=alphabet

    return motifList

def UploadFrommfmd(callback_url, params):
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
        motifList = parse_mfmd_output(params['path'])
        print(motifList)
       
        MSO = {}
        MSO=motifList
        '''MSO['Condition'] = 'Temp'
        MSO['SequenceSet_ref'] = '123'
        MSO['Motifs'] = []
        MSO['Alphabet'] = ['A','C','G','T']
        #MSO['Background'] = MSU.GetBackground()
        #for letter in MSO['Alphabet']:
        #    MSO['Background'][letter] = 0.0
        
        #MSU.parseMotifList(motifList,MSO)'''
        
        '''params['min_len']=22   #put dummy value for min and max len
        params['max_len']=22
        #MSU.CheckLength(motifList,params['min_len'],params['max_len'])
        #MSU.CheckLength(MSO,params['min_len'],params['max_len'])
        
        
        for motif in MSO['Motifs']:
            print()
            for letter in MSO['Alphabet']:
                if len(motif['PWM'][letter]) != len(motif['Iupac_sequence']):
                    print('CAUGHT PWM ERROR HERE')
                    exit(1)
        if 'absolute_locations' in params:
            for motif in MSO['Motifs']:
                for loc in motif['Motif_Locations']:
                    if loc['sequence_id'] in params['absolute_locations']:
                        loc['sequence_id'] = params['contig']
                        absStart = int(params['start'])
                        loc['start'] = absStart
                        loc['end'] = absStart + loc['end']
        print("test2")'''
        
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
        #END UploadFrommfmd

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method UploadFrommfmd return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

