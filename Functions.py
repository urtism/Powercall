import argparse
import subprocess
import json
import os
import textwrap
import random as rm
import datetime
import pysam
import random as r
import regex as re

def prRed(prt): print("\033[91m {}\033[00m" .format(prt))
def prGreen(prt): print("\033[92m {}\033[00m" .format(prt))

def makedirs(dirs):
    for d in dirs:
        if not os.path.exists(d):
			os.makedirs(d)

def panel_check(panel,cfg):
	
	if panel == 'TrusightCardio':
		design = 'Enrichment'
		target_list = cfg['target']['TARGET_TRUSIGHTCARDIO_LIST']
		target_bed = cfg['target']['TARGET_TRUSIGHTCARDIO_BED']
		transcripts_list = cfg['file']['TRANSCR_TRUSIGHTCARDIO']

	elif panel == 'TrusightCancer':
		design = 'Enrichment'
		target_list = cfg['target']['TARGET_TRUSIGHTCANCER_LIST']
		target_bed = cfg['target']['TARGET_TRUSIGHTCANCER_BED']
		transcripts_list = cfg['file']['TRANSCR_TRUSIGHTCANCER']

	elif panel == 'TrusightOne':
		design = 'Enrichment'
		target_list = cfg['target']['TARGET_TRUSIGHTONE_LIST']
		target_bed = cfg['target']['TARGET_TRUSIGHTONE_BED']
		transcripts_list = cfg['file']['TRANSCR_TRUSIGHTONE']

	elif panel == 'BRCAMASTRDx':
		design = 'Amplicon'
		target_list = cfg['target']['TARGET_MULTIPLICOM_BRCA_LIST']
		target_bed = cfg['target']['TARGET_MULTIPLICOM_BRCA_BED']
		transcripts_list = cfg['file']['TRANSCR_BRCA']

	elif panel == 'HTC':
		design = 'Enrichment'
		target_list = cfg['target']['TARGET_HTC_SOPHIA_LIST']
		target_bed = cfg['target']['TARGET_HTC_SOPHIA_BED']
		transcripts_list = cfg['file']['TRANSCR_HTC_SOPHIA']

	elif panel == 'CustomSSQXT':
		design = 'Enrichment'
		target_list = cfg['target']['TARGET_SURESELECT_LIST']
		target_bed = cfg['target']['TARGET_SURESELECT_BED']
		transcripts_list = cfg['file']['TRANSCR_SURESELECT']

	elif panel == 'CustomHPHS':
		design = 'Amplicon'
		target_list = cfg['target']['TARGET_HALOPLEX_LIST']
		target_bed = cfg['target']['TARGET_HALOPLEX_BED']
		transcripts_list = cfg['file']['TRANSCR_HALOPLEX']

	elif panel == 'Custom':
		design = opts.design
		target_list = cfg['target']['TARGET_CUSTOM_LIST']
		target_bed = cfg['target']['TARGET_CUSTOM_BED']
		transcripts_list = cfg['file']['TRANSCR_CUSTOM']

	return design,target_list,target_bed,transcripts_list


def Delete(files):
	for f in files:
		if os.path.isfile(f) or os.path.ispath(f): 
			status = subprocess.call("rm -rf" + f, shell=True)

def Move(files,newdir):
	for f in files:
		if os.path.isfile(f) or os.path.ispath(f): 
			status = subprocess.call("mv " + f + ' ' + newdir, shell=True)