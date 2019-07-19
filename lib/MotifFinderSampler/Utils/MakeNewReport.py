import json
import sys
import os
class MakeNewReport:
  def __init__(self):
      pass
#Pseudo:
#copy report into directory
#write motif set to directory
#pass info back

  def MakeReport(self, htmlDir,motifSet):
      reportPath = '/kb/module/lib/MotifFinderSampler/Utils/Report/*'
      CopyCommand = 'cp -r ' + reportPath + ' ' + htmlDir
      os.system(CopyCommand)
      jsonFName = htmlDir + '/ReportMotif.json'
      with open(jsonFName,'w') as motifjson:
          json.dump(motifSet,motifjson)
      return
