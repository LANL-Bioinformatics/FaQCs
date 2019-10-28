#!/usr/bin/env python3

import sys, os, subprocess, time , re
from xml.dom import minidom
import logging
import log
# 172800s  = 48 hours
QSUB_TIMELIMIT=172800
QSUB_COMMAND=['qsub','-V','-b', 'y', '-j', 'y','-cwd']  
QSTAT_COMMAND=['qstat', '-xml','-j']
QDEL_COMMAND=['qdel', '-j']

ugeLogger = logging.getLogger(__name__)
log.log_init(ugeLogger)
#ugeLogger=log.log_init()
def qsub(cmd, jobname=None, mem='10G', cpu=4, log=None, email=None):
    qsub_cmd = QSUB_COMMAND.copy();

    if jobname is not None:
        qsub_cmd.extend(['-N',jobname])    
    if log is not None:
        qsub_cmd.extend(['-o',log])
    if email is not None:
        qsub_cmd.extend(['-m', 'abe', '-M', email])
    qsub_cmd.extend(['-pe', 'smp',str(cpu)])    
    qsub_cmd.extend(['-l' ,'h_vmem='+ mem + ',mem_free='+mem] )    
        
    qsub_cmd.extend(cmd)
    ugeLogger.info(" ".join(qsub_cmd))
    proc = subprocess.Popen(qsub_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    outs, errs = proc.communicate()
    if proc.returncode != 0:
        ugeLogger.error("Failed %d %s %s" % (proc.returncode, outs, errs))

    m = re.match(r"Your job (\d+)", outs.decode())
    jobid = m.group(1) 

    return jobid

def qwait_all(idlist):
    # take a diction
    while True:
        job_exit=0
        exitIDs=dict()
        for jid in idlist:
            if jid in exitIDs:
                job_exit += 1
            # get job start time if not the job is not exist
            job_start_t, job_submit_t = qstat(jid)
            if job_start_t is None and job_submit_t is None:
                job_exit += 1
                exitIDs[jid] = True
            elif job_start_t is not None and (time.time() - int(job_start_t)) > QSUB_TIMELIMIT:
                job_exit += 1
                qdel(jid)
                ugeLogger.error("Time out on job id %s" % jid)
      
        if len(idlist) == job_exit or not idlist:
            break
        time.sleep(60)
        ugeLogger.debug("Num of finished Job: %d" % job_exit)
    
    return

def qwait(jid):
    
    while True:
        job_exit=0
        job_start_t, job_submit_t = qstat(jid)
        if job_start_t is None and job_submit_t is None:
            job_exit += 1
        elif job_start_t is not None and (time.time() - int(job_start_t)) > QSUB_TIMELIMIT:
            job_exit += 1
            qdel(jid)
            ugeLogger.error("Time out on job id %d" % jid)

        if job_exit == 1 or not jid:
            break
        time.sleep(5)
    return

def qdel(jobid):
    cmd = QDEL_COMMAND.copy()
    cmd.append(jobid)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    outs, errs = proc.communicate()
    if proc.returncode != 0:
        ugeLogger.error("Failed %d %s %s" % (proc.returncode, outs, errs))

def qstat(jobid):
    cmd = QSTAT_COMMAND.copy()
    cmd.append(jobid)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    outs, errs = proc.communicate()
    if proc.returncode != 0:
        ugeLogger.error("Failed %d %s %s" % (proc.returncode, outs, errs))
    
    dom = minidom.parseString(outs.decode())
    start_time_dom=dom.getElementsByTagName('JAT_start_time')
    sub_time_dom=dom.getElementsByTagName('JB_submission_time')
    start_time = start_time_dom[0].firstChild.nodeValue if start_time_dom else None 
    sub_time=sub_time_dom[0].firstChild.nodeValue if sub_time_dom else None
    return start_time, sub_time
