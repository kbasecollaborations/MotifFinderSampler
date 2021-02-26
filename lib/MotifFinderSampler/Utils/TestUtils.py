import os
class TestUtils:
  def __init__(self):
      pass
  def GetGenome(self, targetpath):
      '''

      :param targetpath:
      :return:
      '''
      localgenomepath = '/kb/data/Ptrichocarpa_444_v3.1.fa.assembly.fa'
      command = 'cp ' + localgenomepath + ' ' + targetpath
      os.system(command)
