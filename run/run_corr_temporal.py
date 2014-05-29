#!/usr/bin/env python

import os, sys
import numpy as np
from time import strftime
import timeit
from pandas import read_csv

#basedir         = os.path.abspath("..")
basedir         = "/home2/zarrar/projects/qc"
sys.path.append(os.path.join(basedir, "qclib"))
from qc import cpac_many_subs

# Settings
outdir          = os.path.join(basedir, "corr.qc")
ncores          = 36

# Pipeline/Session
pipeline_dir    = "/data/Projects/CoRR/preproc/output/pipeline_corr_qc_preproc"
info            = read_csv(os.path.join(outdir, "scan_info.csv"), index_col=0)
info.index      = range(info.shape[0])

# Skip the initial runs since they have been run
info            = info.ix[2:,]

for i,row in info.iterrows():
    # i,row = info.iterrows().next()
    session_id  = row.session
    scan_id     = row.scan
    
    print "====="
    print "SESSION: %i" % session_id
    print "SCAN: %i" % scan_id
    print "time: %s" % strftime("%Y-%m-%d %H:%M:%S")
    
    # Subjects
    sublist_file    = "%s/subject_lists/session%02i_scan%02i.txt" % (outdir, session_id, scan_id)
    sublist         = np.loadtxt(sublist_file, dtype=str).tolist()
    if isinstance(sublist, str):
        sublist     = [sublist]
    
    # Check
    outfile = os.path.join(outdir, "qc_temporal_epi_session%02i_scan%02i.csv" % (session_id, scan_id))
    if os.path.exists(outfile):
        out = read_csv(outfile)
        if out.shape[0] == 0:
            print "WARNING: no output exists, rerunning"
        elif out.shape[0] != len(sublist):
            print "ERROR: Number of subjects in output is not the same as input subject list, re-running"
            
            #sublist_file    = "%s/qc_temporal_epi_bad_subs_session%02i_scan%02i.txt" % (outdir, session_id, scan_id)
            #sublist         = np.loadtxt(sublist_file, dtype=str).tolist()
            #if isinstance(sublist, str):
            #    sublist     = [sublist]
            #print "Rerunning %i subjects" % len(sublist)
            
            #continue
        else:
            print "Output already exists, skipping"
            continue    
    
    # exclude 0025770
    #sublist = [ s for s in sublist if s != "0025770" ]
    
    # COMPUTE COMPUTE COMPUTE ...
    start           = timeit.default_timer()
    #(bad_subs, df)  = cpac_many_subs("temporal", sublist[:2], pipeline_dir, 
    #                                  session_id=session_id, scan_id=scan_id, ncores=1)
    (bad_subs, df)  = cpac_many_subs("temporal", sublist, pipeline_dir, 
                                     session_id=session_id, scan_id=scan_id, ncores=ncores)
    stop            = timeit.default_timer()
    duration        = (stop-start)/60.  # in minutes
    print "time elapsed (mins):", duration
    
    # Save Output
    print "...saving output"
    ## df
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outfile = os.path.join(outdir, "qc_temporal_epi_session%02i_scan%02i.csv" % (session_id, scan_id))
    df.to_csv(outfile)
    ## bad subjects
    outfile = os.path.join(outdir, "qc_temporal_epi_bad_subs_session%02i_scan%02i.txt" % (session_id, scan_id))
    np.savetxt(outfile, bad_subs, fmt="%s")
    
    print "====="
