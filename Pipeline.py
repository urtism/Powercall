import subprocess
import alignment
import preprocessing


def Pipeline_Germline_Multisample(workflow,samplesheet,design,panel,dirs,cfg,opts,target_list,target_bed,transcripts_list):

	status = subprocess.call("cat "+ dirs['logo']+'/logo_cmg.txt', shell=True)

	if 'A' in workflow:

		samples = alignment.ReadSampleSheet(samplesheet,'Germline')
		alignment_log = open(dirs['log'] + '/Alignment.log','w+')
		samplesheet = dirs['log'] + '/preprocessing.samplesheet'
		new_samplesheet = open(samplesheet,'w+')

		for sample in samples.keys():
			if panel == 'CustomHPHS':
				sample_name,fq1,fq2,fqI2 = samples[sample][0]
				bam = alignment.Alignment(panel,sample_name,fq1,fq2,fqI2,dirs,cfg,opts,alignment_log)
			else:
				sample_name,fq1,fq2 = samplesheet[sample][0]
				bam = alignment.Alignment(panel,sample_name,fq1,fq2,None,dirs,cfg,opts,alignment_log)

			new_samplesheet.write(sample_name+'\t'+bam+'\n')
		alignment_log.close()
		new_samplesheet.close()

	if 'R' in workflow or 'M' in workflow or 'I' in workflow or 'B' in workflow:

		samples = preprocessing.ReadSampleSheet(samplesheet,'Germline')
		preprocessing_log = open(dirs['log'] + '/Preprocessing.log','w+')
		samplesheet = dirs['log'] + '/variantcalling.samplesheet'
		new_samplesheet = open(samplesheet,'w+')

		for sample in samples.keys():
			
			sample_name,bam = samplesheet[sample][0]
			bam = preprocessing.Preprocessing(panel,sample_name,fq1,fq2,None,dirs,cfg,opts,alignment_log)

			new_samplesheet.write(sample_name+'\t'+bam+'\n')

		preprocessing_log.close()
		new_samplesheet.close()
	
	if 'V' in workflow:

		VARIANT_CALLING_GERMLINE $CFG

	if 'F' in workflow:
		Features_extraction_germline $CFG

	if 'E' in workflow:




def pipeline.Pipeline_Germline_Singlesample(workflow,samplesheet,design,panel,dirs,cfg,opts,target_list,target_bed,transcripts_list):

	status = subprocess.call("cat "+ dirs['logo']+'/logo_cmg.txt', shell=True)

	if 'A' in workflow:
		samples = alignment.ReadSampleSheet(samplesheet,'Germline')
		ALLINEAMENTO $CFG

	if 'R' in workflow or 'M' in workflow or 'I' in workflow or 'B' in workflow:

		PREPROCESSING $CFG
	
	if 'V' in workflow:

		VARIANT_CALLING_GERMLINE $CFG

	if 'F' in workflow:
		Features_extraction_germline $CFG

	if 'E' in workflow:


def pipeline.Pipeline_Somatic_Case_Control(workflow,samplesheet,design,panel,dirs,cfg,opts,target_list,target_bed,transcripts_list):

	status = subprocess.call("cat "+ dirs['logo']+'/logo_cmg.txt', shell=True)

	if 'A' in workflow:
		samples = alignment.ReadSampleSheet(samplesheet,'Somatic_Case_Control')
		ALLINEAMENTO $CFG

	if 'R' in workflow or 'M' in workflow or 'I' in workflow or 'B' in workflow:

		PREPROCESSING $CFG
	
	if 'V' in workflow:

		VARIANT_CALLING_GERMLINE $CFG

	if 'F' in workflow:
		Features_extraction_germline $CFG

	if 'E' in workflow:

def pipeline.Pipeline_Somatic(workflow,samplesheet,design,panel,dirs,cfg,opts,target_list,target_bed,transcripts_list):

	status = subprocess.call("cat "+ dirs['logo']+'/logo_cmg.txt', shell=True)

	if 'A' in workflow:
		samples = alignment.ReadSampleSheet(samplesheet,'Somatic')
		ALLINEAMENTO $CFG

	if 'R' in workflow or 'M' in workflow or 'I' in workflow or 'B' in workflow:

		PREPROCESSING $CFG
	
	if 'V' in workflow:

		VARIANT_CALLING_GERMLINE $CFG

	if 'F' in workflow:
		Features_extraction_germline $CFG

	if 'E' in workflow: