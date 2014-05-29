#!/usr/bin/env python

import os
import numpy as np
from time import strftime
import timeit

#basedir         = os.path.abspath("..")
basedir         = "/home2/zarrar/projects/qc"
sys.path.append(os.path.join(basedir, "qclib"))
from qc import cpac_many_subs

# Settings
outdir          = "corr.qc"
ncores          = 20

# Pipeline/Session
session_ids     = range(5,11) # ran 1-4
scan_id         = 1
pipeline_dir    = "/data/Projects/CoRR/preproc/output/pipeline_corr_qc_preproc"

for session_id in session_ids:
    print "====="
    print "SESSION: %i" % session_id
    print "time: %s" % strftime("%Y-%m-%d %H:%M:%S")
    
    # Subjects
    sublist_file    = "%s/subject_lists/session%02i_scan%02i.txt" % (outdir, session_id, scan_id)
    sublist         = np.loadtxt(sublist_file, dtype=str).tolist()
    
    # COMPUTE COMPUTE COMPUTE ...
    start           = timeit.default_timer()
    #(bad_subs, df)  = cpac_many_subs("spatial", sublist[:10], pipeline_dir, session_id=session_id, ncores=1)
    #(bad_subs, df)  = cpac_many_subs("spatial", ["0026048"], pipeline_dir, session_id=session_id, ncores=1)
    (bad_subs, df)  = cpac_many_subs("spatial", sublist, pipeline_dir, session_id=session_id, ncores=ncores)
    stop            = timeit.default_timer()
    duration        = (stop-start)/60.  # in minutes
    print "time elapsed (mins):", duration
    
    
    # Save Output
    print "...saving output"
    ## df
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    outfile = os.path.join(outdir, "qc_spatial_session%02i.csv" % session_id)
    df.to_csv(outfile)
    ## bad subjects
    outfile = os.path.join(outdir, "qc_spatial_bad_subs_session%02i.txt" % session_id)
    np.savetxt(outfile, bad_subs, fmt="%s")
    
    print "====="
