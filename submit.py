#!/cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw/CMSSW_13_3_3/external/el8_amd64_gcc12/bin/python3

import sys
import subprocess
import datetime
import numpy as np
import glob

def checkRunning(processes) :
  nRunning = 0
  for proc, name in processes :
    status = proc.poll()
    if status==None : nRunning+=1
  return nRunning
    
dataDir='/scratch_cmsse/konec/data/2023D_ParkingDoubleMuonLowMass/'
lsCommand='ls -1 '+dataDir+'|grep root'
print ('command: ',lsCommand)
dir=subprocess.Popen(lsCommand, stdout=subprocess.PIPE,shell=True,text=True)
lsOutput=dir.communicate()[0]
files=[]
for f in lsOutput.split():
  print(dataDir+f)
  files.append(dataDir+f)
print ('number of files: ',len(files))

nJobs = 216
maxRunningJobs = 8 
files_splitted = np.array_split(files, nJobs)
print ('number of files: ',len(files),', submitting in', nJobs,' jobs, max running jobs: ',maxRunningJobs)

myProc=[]
for fs in range(nJobs):
  while (checkRunning(myProc) >= maxRunningJobs) : subprocess.run(['sleep','5'])

  print ('nRunning is:', checkRunning(myProc),' submitting job number:',fs,'\n', files_splitted[fs])
  jobId =str(fs).zfill(3)
  execCommand = ['cmsRun','../analysis_cwKSjob.py', jobId, str(files_splitted[fs]).strip('[]')]
  print ('execCommand #',fs,' is: ', execCommand)
  p=subprocess.Popen(execCommand,stdout=open('out_'+jobId+'.txt','w'), stderr=subprocess.STDOUT)
  subprocess.run(['sleep','5'])
  myProc.append( (p, jobId) )

print ('...all submitted.')
while (checkRunning(myProc)) : subprocess.run(['sleep','30'])
print ('all finished:')
for proc,name in myProc :
  status=proc.poll()
  print ('status of', proc, name, 'is: ',status) 
  subprocess.run(['sleep','1'])

matching_files = glob.glob('histos_*.root')
print( matching_files )
haddCommand = ['hadd', '-f', 'for5_27_05.root'] + matching_files
subprocess.run(haddCommand)

sys.exit(0)
