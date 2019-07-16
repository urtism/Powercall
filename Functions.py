import argparse
import subprocess
from ruffus import *
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


def init_dirs(work_dir,opts):

	dirs = dict()
	dirs['work'] = work_dir
	dirs['storage'] = work_dir + '/STORAGE/'+opts.run_id
	dirs['out'] = work_dir + '/OUTPUT/'+opts.run_id
	dirs['log'] = work_dir + '/LOGS'
	dirs['delete'] = work_dir + '/DELETE'
	dirs['alignment'] = work_dir + '/ALIGNMENT'
	dirs['preprocessing'] = work_dir + '/PREPROCESSING'
	dirs['variantcalling'] = work_dir + '/VARIANTCALLING'
	dirs['gvcf'] = work_dir + '/VARIANTCALLING/GVCF'
	dirs['featsextract'] = work_dir + '/FEATURES_EXTRACTION'
	dirs['annotation'] = work_dir + '/ANNOTATION'
	dirs['script'] = os.path.dirname(os.path.abspath(__file__)) + '/scripts/'
	dirs['logo'] = os.path.dirname(os.path.abspath(__file__)) + '/LOGOS'
	
	return dirs

def makedirs(dirs):
	for d in dirs:
		if not os.path.exists(d):
			os.makedirs(d)

def panel_check(panel,cfg):
	
	if panel == 'TrusightCardio':
		design = 'Enrichment'
		target_list = cfg['target']['TARGET_TRUSIGHTCARDIO_LIST']
		target_bed = cfg['target']['TARGET_TRUSIGHTCARDIO_BED']
		transcripts_list = cfg['files']['TRANSCR_TRUSIGHTCARDIO']

	elif panel == 'TrusightCancer':
		design = 'Enrichment'
		target_list = cfg['target']['TARGET_TRUSIGHTCANCER_LIST']
		target_bed = cfg['target']['TARGET_TRUSIGHTCANCER_BED']
		transcripts_list = cfg['files']['TRANSCR_TRUSIGHTCANCER']

	elif panel == 'TrusightOne':
		design = 'Enrichment'
		target_list = cfg['target']['TARGET_TRUSIGHTONE_LIST']
		target_bed = cfg['target']['TARGET_TRUSIGHTONE_BED']
		transcripts_list = cfg['files']['TRANSCR_TRUSIGHTONE']

	elif panel == 'BRCAMASTRDx':
		design = 'Amplicon'
		target_list = cfg['target']['TARGET_MULTIPLICOM_BRCA_LIST']
		target_bed = cfg['target']['TARGET_MULTIPLICOM_BRCA_BED']
		transcripts_list = cfg['files']['TRANSCR_BRCA']

	elif panel == 'HTC':
		design = 'Enrichment'
		target_list = cfg['target']['TARGET_HTC_SOPHIA_LIST']
		target_bed = cfg['target']['TARGET_HTC_SOPHIA_BED']
		transcripts_list = cfg['files']['TRANSCR_HTC_SOPHIA']

	elif panel == 'CustomSSQXT':
		design = 'Enrichment'
		target_list = cfg['target']['TARGET_SURESELECT_LIST']
		target_bed = cfg['target']['TARGET_SURESELECT_BED']
		transcripts_list = cfg['files']['TRANSCR_SURESELECT']

	elif panel == 'CustomHPHS':
		design = 'Amplicon'
		target_list = cfg['target']['TARGET_HALOPLEX_LIST']
		target_bed = cfg['target']['TARGET_HALOPLEX_BED']
		transcripts_list = cfg['files']['TRANSCR_HALOPLEX']

	elif panel == 'Custom':
		design = opts.design
		target_list = cfg['target']['TARGET_CUSTOM_LIST']
		target_bed = cfg['target']['TARGET_CUSTOM_BED']
		transcripts_list = cfg['files']['TRANSCR_CUSTOM']

	return design,target_list,target_bed,transcripts_list

def ReadSampleSheet(samplesheet,analysis,panel,step):
	ssheet = dict()
	if step == 'Alignment':
		with open(samplesheet,'r') as ss:
			for line in ss:
				line=line.rstrip()
				if panel == 'CustomHPHS':
					if analysis == 'Germline' or analysis == 'Somatic':
						sample_name = line.split('\t')[0]
						fq1 = line.split('\t')[1]
						fq2 = line.split('\t')[2]
						fqI2 = line.split('\t')[3]
						ssheet[sample_name] = [[sample_name,fq1,fq2,fqI2]]

					elif analysis == 'Somatic_Case_Control':
						gsample_name = line.split('\t')[0]
						gfq1 = line.split('\t')[1]
						gfq2 = line.split('\t')[2]
						gfqI2 = line.split('\t')[3]
						ssample_name = line.split('\t')[4]
						sfq1 = line.split('\t')[5]
						sfq2 = line.split('\t')[6]
						sfqI2 = line.split('\t')[7]
						ssheet[ssample_name] = [[gsample_name,gfq1,gfq2,gfqI2],[ssample_name,sfq1,sfq2,sfqI2]]
				else:
					if analysis == 'Germline' or analysis == 'Somatic':
						sample_name = line.split('\t')[0]
						fq1 = line.split('\t')[1]
						fq2 = line.split('\t')[2]
						ssheet[sample_name] = [[sample_name,fq1,fq2]]

					elif analysis == 'Somatic_Case_Control':
						gsample_name = line.split('\t')[0]
						gfq1 = line.split('\t')[1]
						gfq2 = line.split('\t')[2]
						ssample_name = line.split('\t')[3]
						sfq1 = line.split('\t')[4]
						sfq2 = line.split('\t')[5]
						ssheet[ssample_name] = [[gsample_name,gfq1,gfq2],[ssample_name,sfq1,sfq2]]
	
	elif step == 'Preprocessing':
		with open(samplesheet,'r') as ss:
			for line in ss:
				line=line.rstrip()
				if analysis == 'Germline' or analysis == 'Somatic':
					sample_name = line.split('\t')[0]
					bam = line.split('\t')[1]
					ssheet[sample_name] = [[sample_name,bam]]

				elif analysis == 'Somatic_Case_Control':
					gsample_name = line.split('\t')[0]
					gbam = line.split('\t')[1]
					ssample_name = line.split('\t')[2]
					sbam = line.split('\t')[3]
					ssheet[ssample_name] = [[gsample_name,gbam],[ssample_name,sbam]]

	elif step == 'VariantCalling':
		with open(samplesheet,'r') as ss:
			for line in ss:
				line=line.rstrip()
				if analysis == 'Germline' or analysis == 'Somatic':
					sample_name = line.split('\t')[0]
					bam = line.split('\t')[1]
					ssheet[sample_name] = [[sample_name,bam]]

				elif analysis == 'Somatic_Case_Control':
					gsample_name = line.split('\t')[0]
					gbam = line.split('\t')[1]
					ssample_name = line.split('\t')[2]
					sbam = line.split('\t')[3]
					ssheet[ssample_name] = [[gsample_name,gbam],[ssample_name,sbam]]

	elif step == 'Featuresextraction':
		with open(samplesheet,'r') as ss:
			for line in ss:
				line=line.rstrip()
				if analysis == 'Germline':
					sample_name = line.split('\t')[0]
					try:
						vcf_gatk = line.split('\t')[1]
						vcf_freebayes = line.split('\t')[2]
						vcf_varscan = line.split('\t')[3]
						ssheet[sample_name] = [[sample_name,vcf_gatk,vcf_freebayes,vcf_varscan]]
					except:
						vcf_merge = line.split('\t')[1]
						ssheet[sample_name] = [[sample_name,vcf_merge]]
				elif analysis == 'Somatic_Case_Control':
					gsample_name = line.split('\t')[0]
					ssample_name = line.split('\t')[1]
					mutect_vcf,vardict_vcf,varscan_vcf = line.split('\t')[2:]
					ssheet[ssample_name] = [[gsample_name,ssample_name,mutect_vcf,vardict_vcf,varscan_vcf]]
				elif analysis == 'Somatic':
					ssample_name = line.split('\t')[0]
					mutect_vcf = line.split('\t')[1]
					ssheet[ssample_name] = [[ssample_name,mutect_vcf]]

	elif step == 'Annotation':
		with open(samplesheet,'r') as ss:
			for line in ss:
				line=line.rstrip()
				if analysis == 'Germline':
					name = line.split('\t')[0]
					vcf = line.split('\t')[1]
					samples = line.split('\t')[2]
					ssheet[name] = [[name,vcf,samples]]
					#ssheet[name] = [[name,vcf]]
				elif analysis == 'Somatic_Case_Control' or analysis == 'Somatic':
					name = line.split('\t')[0]
					vcf = line.split('\t')[1]
					tsv = line.split('\t')[2]
					ssheet[name] = [[name,vcf,tsv]]
			
	return ssheet


def Delete(files):
	for f in files:
		if os.path.isfile(f): 
			status = subprocess.call("rm -rf" + f, shell=True)
		elif os.path.ispath(f):
			status = subprocess.call("rm -rf" + f, shell=True)

def Move(files,newdir):
	for f in files:
		if os.path.isfile(f): 
			status = subprocess.call("mv " + f + ' ' + newdir, shell=True)
		elif os.path.ispath(f):
			status = subprocess.call("mv " + f + ' ' + newdir, shell=True)

def Copy(files,newdir):
	for f in files:
		if os.path.isfile(f): 
			status = subprocess.call("cp " + f + ' ' + newdir, shell=True)
		else:
			print files