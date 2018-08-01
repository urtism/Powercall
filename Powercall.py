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
import Pipeline as pipeline
import Functions as f


if __name__ == '__main__':

	parser = argparse.ArgumentParser('This is the most cool Pipeline you have ever seen')
	parser.add_argument('-v', '--version', action='version', version='Powercall v3.2.0')

	parser.add_argument('-c', '--cfg', help="Configuration file in json format", required=True)
	parser.add_argument('-s', '--samplesheet', help="Samplesheet that contains filepaths and samples associated", required=True)
	parser.add_argument('-d', '--design', help="Design of the NGS experiment [Enrichment,Amplicon].",choices=['Enrichment','Amplicon'])
	parser.add_argument('-p', '--panel', help="Type of Panel used for the analysis [TrusightCardio (illumina), TrusightCancer (illumina), TrusightOne (illumina), BRCAMASTRDx (MULTIPLICOM), HTC (SophiaGenetics), CustomSSQXT (Custom SureSelect QXT Agilent), CustomHPHS (Custom HaloPlex HS Agilent), Custom (Other Custom panel)]",
		choices="[TrusightCardio, TrusightCancer, TrusightOne, BRCAMASTRDx, HTC, CustomSSQXT, CustomHPHS, Custom]")
	parser.add_argument('-a', '--analysis', help="Type of analysis [Germline (Multisample), GermlineSS (Singlesample), SomaticCC (Case-Control), Somatic (only Case)].", choices=['Germline', 'GermlineSS','SomaticCC','Somatic'], default='Germline')
	parser.add_argument('-id','--run_id', help="Run id [YYYYMMDD_Run_NRUN_PANEL}", required=True)
	parser.add_argument('-w', '--workdir', help="Working Directory. Use this option if you want to work in a precise directory. Default: '~/NGS_ANALYSIS/run_id'")
	parser.add_argument('--workflow', help="String that indicates which steps the analysis must do: A-> Alignment, R-> AddOrReplaceReadGroups, M-> MarkDuplicates, I-> IndelRealigner, B-> BaseRecalibrator, V-> Variant Calling, F-> Features Extraction, E-> Annotation. Es: --start AMIBVFE indicates all steps (like --start ALL), --start MIBV starts from MarkDuplicates and ends to Variant Calling. Default: ALL ", default='ALL')

		
	start_time = datetime.datetime.now()
	
	global opts
	opts = parser.parse_args()

	dirs=dict()

	if opts.cfg != None:
		cfg = json.loads((open(opts.cfg).read()).encode('utf8'))
	else:
		cfg = json.loads((open(os.path.dirname(os.path.abspath(__file__)) + '/CFG/Powercall.default.cfg.json').read()).encode('utf8'))

	script_dir = os.path.dirname(os.path.abspath(__file__)) + '/scripts/'
	if opts.workdir != None:
		work_dir = opts.workdir
	else:
		work_dir = '~/NGS_ANALYSIS/' + opts.run_id

	logos_dir = os.path.dirname(os.path.abspath(__file__)) + '/LOGOS'
	storage_dir = work_dir + '/STORAGE/'+opts.run_id
	out_dir = work_dir + '/OUTPUT/'+opts.run_id
	delete_dir = work_dir + '/DELETE'
	log_dir = work_dir + '/LOGS'
	align_dir = work_dir + '/ALIGNMENT'
	preprocessing_dir = work_dir + '/PREPROCESSING'
	variantcalling_dir = work_dir + '/VARIANTCALLING'
	gcvf_dir = work_dir + '/VARIANTCALLING/GVCF'
	featuresextraction_dir = work_dir + '/FEATURES_EXTRACTION'
	annotation_dir = work_dir + '/ANNOTATION'

	dirs['work'] = work_dir
	dirs['storage'] = storage_dir
	dirs['out'] = out_dir
	dirs['log'] = log_dir
	dirs['delete'] = delete_dir
	dirs['logo'] = logos_dir
	dirs['alignment'] = align_dir
	dirs['preprocessing'] = preprocessing_dir
	dirs['variantcalling'] = variantcalling_dir
	dirs['gvcf'] = gcvf_dir
	dirs['featsextract'] = featuresextraction_dir
	dirs['annotation'] = annotation_dir
	dirs['script'] = script_dir
	
	f.makedirs([work_dir, work_dir+'/STORAGE', work_dir+'/OUTPUT', storage_dir, out_dir, delete_dir,log_dir])

	design,target_list,target_bed,transcripts_list = f.panel_check(opts.panel,cfg)

	if opts.design != None:
		design = opts.design

	workflow = opts.workflow
	if workflow == 'ALL':
		workflow = 'ARMIBVFE'

	if design == 'Amplicon':
		workflow=re.sub('MIB','',workflow)

	if opts.analysis == 'Germline':
		pipeline.Pipeline_Germline_Multisample(workflow,opts.samplesheet,design,opts.panel,dirs,cfg,opts,target_list,target_bed,transcripts_list)

	if opts.analysis == 'GermlineSS':
		pipeline.Pipeline_Germline_Singlesample(workflow,opts.samplesheet,design,opts.panel,dirs,cfg,opts,target_list,target_bed,transcripts_list)

	elif opts.analysis == 'SomaticCC':
		pipeline.Pipeline_Somatic_Case_Control(workflow,opts.samplesheet,design,opts.panel,dirs,cfg,opts,target_list,target_bed,transcripts_list)

	elif opts.analysis == 'Somatic':
		pipeline.Pipeline_Somatic(workflow,opts.samplesheet,design,dirs,opts.panel,cfg,opts,target_list,target_bed,transcripts_list)


	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)
	print "\nTime elapsed: %d min, %d sec" % elapsed_time