###
# Import various things
###

import sys
basedir         = "/home2/zarrar/projects/qc"
sys.path.append(os.path.join(basedir, "qclib"))

import os
import numpy as np
from time import strftime
import timeit
from pandas import read_csv

from qc import cpac_many_subs
from qc import cpac_qc_temporal

from qc import cpac_qc_spatial
from qc import CpacPaths2
from spatial_qc import efc
from qc import load_image
from tempfile import mkstemp


###
# Data Analysis
###

# Replicate error
tmp = cpac_qc_temporal("0025770", "/data/Projects/CoRR/preproc/output/pipeline_corr_qc_preproc", session_id=2, scan_id=1)        

# Useful inputs
pipeline_dir="/data/Projects/CoRR/preproc/output/pipeline_corr_qc_preproc"
subject_id="0026019" 
session_id=1
scan_id=1 

# Get paths
cpac_paths  = CpacPaths2(pipeline_dir, subject_id, session_id, None)
fg_path     = cpac_paths.extract("functional_brain_mask", "_scan_rest_%i*_rest/*.nii.gz" % scan_id)
# Create mean EPI with background noise
func_path   = cpac_paths.extract("motion_correct", "_scan_rest_%i*_rest/*.nii.gz" % scan_id)
_,anat_path = mkstemp(suffix=".nii.gz", prefix="tmp_mean_epi_")
cmd         = "fslmaths %s -Tmean %s" % (func_path, anat_path)
os.system(cmd)

# Read in data
anat_data       = load_image(anat_path)

# Check the scale and data type (this might be an issue)
import nibabel as nib
img             = nib.load(anat_path)
hdr             = img.get_header()
hdr.get_slope_inter()   # (1.0, 0.0)
img.get_data_dtype()    # float32

# this gives a nan
res = efc(anat_data)

# let's try it manually
# Calculate the maximum value of the EFC (which occurs any time all voxels have the same value)
efc_max = 1.0 * np.prod(anat_data.shape) * (1.0 / np.sqrt(np.prod(anat_data.shape))) * \
            np.log(1.0 / np.sqrt(np.prod(anat_data.shape)))
sh
# Calculate the total image energy
b_max   = np.sqrt((anat_data**2).sum())

# Calculate EFC (add 1e-16 to the image data to keep log happy)
efc2     = (1.0 / efc_max) * np.sum((anat_data / b_max) * np.log((anat_data + 1e-16) / b_max))
## will get nan


###
# Getting to the bottom of this
###

# Will see that values that cause nan are all negative (note: anat_data is actually mean_func)
anat_data[np.isnan(np.log((anat_data+1e-16)/b_max))]

# Examine the functional data
img = nib.load(func_path)
dat = img.get_data()
hdr = img.get_header()
hdr.get_slope_inter()   # no scale
hdr.get_data_offset()
hdr.get_data_dtype()    # int

# Look into the motion corrected functional data for points causing nan
# Negative too
dat[np.isnan(np.log((anat_data+1e-16)/b_max))].mean(1)  # show mean across time
dat[np.isnan(np.log((anat_data+1e-16)/b_max))][0,:]     # visualize some offending point time-series
dat[np.isnan(np.log((anat_data+1e-16)/b_max))][1,:]
dat[np.isnan(np.log((anat_data+1e-16)/b_max))][2,:]
dat[np.isnan(np.log((anat_data+1e-16)/b_max))][4,:]

# Get the raw data as well and look at it
raw_path = "/home/data/Incoming/CORR/ExtractedArchives/Organized_Data/UTAH_org/0026019/session_1/rest_1/rest.nii.gz"
raw_dat = load_image(raw_path)
raw_dat.dtype   # also int
raw_dat[np.isnan(np.log((anat_data+1e-16)/b_max))].mean(1)  # not all are negative

# 

# Number of voxels with negative time-points
tmp0 = (raw_dat<0).any(3)
tmp1 = (dat<0).any(3)
tmp0.sum()
tmp1.sum()

# Visualize offending time-series
raw_dat[tmp0]
raw_dat[tmp0][0]
raw_dat[tmp0][1]
