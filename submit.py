#!/usr/bin/python


import sys
import os
from commands import getstatusoutput
import datetime


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


def prepareJob (jobID):
    filename = 'job_' + jobID
    f = open (filename, 'w')
    f.write ('cd /afs/cern.ch/user/g/govoni/w0/shashlik/Shashlik\n')
    f.write ('source setup_lxplus.sh\n')
    f.write ('/afs/cern.ch/user/g/govoni/geant4_workdir/bin/Linux-g++/Shashlik config.cfg out_' + jobID + '\n')
    f.close ()
    getstatusoutput ('chmod 755 job_' + jobID)
    getstatusoutput ('bsub -q 8nh job_' + jobID)


# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


if __name__ == '__main__':

    if len (sys.argv) < 2 : 
        print 'usage :', sys.argv[0], 'jobsNumber','\n'
        exit (1)
    jobsNum = int (sys.argv[1])
        
    print 'submitting', jobsNum, 'jobs'
        
    for i in range (0, jobsNum):
        prepareJob (str (i))    

