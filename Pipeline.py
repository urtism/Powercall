import subprocess
import tools
import datetime
import Functions as f


def Alignment(panel,sample_name,fq1,fq2,fqI2,dirs,cfg,opts,log,workdir):
	start_time = datetime.datetime.now()
	print 'Sample: '+sample_name
	f.makedirs([dirs['alignment']])

	workdir = dirs['alignment']

	if panel == 'CustomSSQXT':

		fq1,fq2 = tools.SureCallTrimmer(cfg['tools']['SURECALLTRIMMER']['path'],cfg['tools']['SURECALLTRIMMER']['ram'],fq1,fq2,'-hs',dirs['workdir']+'/TRIM_FASTQ',log)
		sam = tools.Bwa_mem(cfg['tools']['BWA']['path'],cfg['tools']['BWA']['threads'],sample_name,fq1,fq2,cfg['reference']['FASTA'],log,workdir)
		bam = tools.SamFormatConverter(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sam,log,workdir)
		mbbam = tools.LocatIt(cfg['tools']['LOCATIT']['path'],cfg['tools']['LOCATIT']['ram'],bam,fqI2,log,workdir)
		sort = tools.SortSam(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],mbbam,log,workdir)
		tools.BuildBamIndex(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sort,log)

		#f.Delete([sam,bam,mbbam])
		f.Move([sam,bam,mbbam],dirs['delete'])
	elif panel == 'CustomHPHS':

		fq1,fq2 = tools.SureCallTrimmer(cfg['tools']['SURECALLTRIMMER']['path'],cfg['tools']['SURECALLTRIMMER']['ram'],fq1,fq2,'-qxt',dirs['workdir']+'/TRIM_FASTQ',log,)
		sam = tools.Bwa_mem(cfg['tools']['BWA']['path'],cfg['tools']['BWA']['threads'],sample_name,fq1,fq2,cfg['reference']['FASTA'],log,workdir)
		bam = tools.SamFormatConverter(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sam,log,workdir)
		sort = tools.SortSam(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],bam,log,workdir)
		tools.BuildBamIndex(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sort,log)

		#f.Delete([sam,bam])
		f.Move([sam,bam],dirs['delete'])
	else:

		sam = tools.Bwa_mem(cfg['tools']['BWA']['path'],cfg['tools']['BWA']['threads'],sample_name,fq1,fq2,cfg['reference']['FASTA'],log,workdir)
		bam = tools.SamFormatConverter(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sam,log,workdir)
		sort = tools.SortSam(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],bam,log,workdir)
		tools.BuildBamIndex(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sort,log)

		#f.Delete([sam,bam])
		f.Move([sam,bam],dirs['delete'])
	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)
	print 'Sample: '+sample_name+ ' -> Done in: %d min, %d sec' % elapsed_time
	return sort


def Preprocessing(panel,workflow,target_list,sample_name,bam,dirs,cfg,opts,log,workdir):

	print 'Sample: ' + sample_name
	start_time = datetime.datetime.now()

	if 'R' in workflow:

		ARRG_bam = tools.AddOrReplaceReadGroups(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],bam,sample_name,panel,opts.run_id,log,workdir)
		tools.BuildBamIndex(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],ARRG_bam,log)
		bam = ARRG_bam

	if 'M' in workflow:

		MD_bam = tools.MarkDuplicates(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],bam,log,workdir)
		tools.BuildBamIndex(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],MD_bam,log)
		#f.Delete([sam,bam])
		f.Move([bam],dirs['delete'])
		bam = MD_bam

	if 'I' in workflow:

		IR_bam = tools.IndelRealigner(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],bam,cfg['reference']['FASTA'],cfg['database']['MILLS'],target_list,log,workdir)
		tools.BuildBamIndex(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],IR_bam,log)
		#f.Delete([sam,bam])
		f.Move([bam],dirs['delete'])
		bam = IR_bam

	if 'B' in workflow:

		BR_bam = tools.BaseRecalibrator(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],bam,cfg['reference']['FASTA'],cfg['database']['DBSNP'],cfg['database']['MILLS'],target_list,log,workdir)
		#f.Delete([sam,bam])
		f.Move([bam],dirs['delete'])
		bam = BR_bam

	tools.BuildBamIndex((cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],bam,log))

	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)
	print 'Sample: '+sample_name+ ' -> Done in: %d min, %d sec' % elapsed_time
	return bam
	#f.Delete([sam,bam])


def Pipeline_Germline_Multisample(workflow,samplesheet,design,panel,dirs,cfg,opts,target_list,target_bed,transcripts_list):

	status = subprocess.call("cat "+ dirs['logo']+'/logo_cmg.txt', shell=True)

	if 'A' in workflow:

		print "ALIGNMENT"
		samples = f.ReadSampleSheet(samplesheet,'Germline',panel,'Alignment')
		alignment_log = open(dirs['log'] + '/Alignment.log','w+')
		samplesheet = dirs['log'] + '/preprocessing.samplesheet'
		new_samplesheet = open(samplesheet,'w+')

		f.makedirs([dirs['alignment']])
		workdir = dirs['alignment']

		for sample in samples.keys():

			if panel == 'CustomHPHS':
				sample_name,fq1,fq2,fqI2 = samples[sample][0]
				bam = Alignment(panel,sample_name,fq1,fq2,fqI2,dirs,cfg,opts,alignment_log,workdir)
			else:
				sample_name,fq1,fq2 = samples[sample][0]
				bam = Alignment(panel,sample_name,fq1,fq2,None,dirs,cfg,opts,alignment_log,workdir)

			new_samplesheet.write(sample_name+'\t'+bam+'\n')
		alignment_log.close()
		new_samplesheet.close()

	if 'R' in workflow or 'M' in workflow or 'I' in workflow or 'B' in workflow:
		print "\nPREPROCESSING"
		samples = f.ReadSampleSheet(samplesheet,'Germline',panel,'Preprocessing')
		preprocessing_log = open(dirs['log'] + '/Preprocessing.log','w+')
		samplesheet = dirs['log'] + '/variantcalling.samplesheet'
		new_samplesheet = open(samplesheet,'w+')

		f.makedirs([dirs['preprocessing']])
		workdir = dirs['preprocessing']

		for sample in samples.keys():
			sample_name,bam = samples[sample][0]
			bam = Preprocessing(panel,workflow,target_list,sample_name,bam,dirs,cfg,opts,preprocessing_log,workdir)

			new_samplesheet.write(sample_name+'\t'+bam+'\n')

		preprocessing_log.close()
		new_samplesheet.close()
	
	if 'V' in workflow:
		print "VARIANT CALLING"

		samples = f.ReadSampleSheet(samplesheet,'Germline',panel,'Preprocessing')
		variantcalling_log = open(dirs['log'] + '/VariantCalling.log','w+')
		samplesheet = dirs['log'] + '/Featuresextraction.samplesheet'
		
		new_samplesheet = open(samplesheet,'w+')

		workdir = dirs['variantcalling']
		f.makedirs([dirs['variantcalling'],workdir+'/GVCF'])

		bam_list = workdir+'/bams.list'
		sample_list = workdir+'/samples.list'
		gvcf_list = workdir+'/gvcf.list'

		bam_l = open(bam_list,'w+')
		sample_l = open(sample_list,'w+')
		gvcf_l = open(gvcf_list,'w+')

		for sample in samples.keys():
			sample_name,bam = samples[sample][0]

			gvcf = tools.HaplotypeCaller(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],bam,sample_name,cfg['reference']['FASTA'],target_list,variantcalling_log,workdir+'/GVCF')

			bam_l.write(bam+'\n')
			sample_l.write(sample_name+'\n')
			gvcf_l.write(gvcf+'\n')

		bam_l.close()
		sample_l.close()
		gvcf_l.close()

		name = opts.run_id

		gatk_vcf = tools.GenotypeGVCFs(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],name,gvcf_list,cfg['reference']['FASTA'],variantcalling_log,workdir)
		gatk_vcf = tools.header_fix(dirs['script']+'merge_vcf.py',gatk_vcf,'G',variantcalling_log)
		gatk_vcf = tools.vcf_norm(cfg['tools']['BEDTOOLS']['path'],gatk_vcf,cfg['reference']['FASTA'],variantcalling_log)
		
		freeb_vcf = tools.FreeBayes(cfg['variantcaller']['FREEBAYES']['path'],name,None,bam_list,cfg['reference']['FASTA'],target_bed,variantcalling_log,workdir)
		freeb_vcf = tools.header_fix(dirs['script']+'merge_vcf.py',freeb_vcf,'F',variantcalling_log)
		freeb_vcf = tools.vcf_norm(cfg['tools']['BEDTOOLS']['path'],freeb_vcf,cfg['reference']['FASTA'],variantcalling_log)
		
		mpileup = tools.mpileup(name,cfg['reference']['FASTA'],None,None,bam_list,target_list,variantcalling_log,workdir)
		varscan_snp_vcf = tools.VarScan_mpileup2snp(cfg['variantcaller']['VARSCAN']['path'],cfg['variantcaller']['VARSCAN']['ram'],mpileup,sample_list,cfg['reference']['FASTA'],variantcalling_log,workdir)
		varscan_indel_vcf = tools.VarScan_mpileup2indel(cfg['variantcaller']['VARSCAN']['path'],cfg['variantcaller']['VARSCAN']['ram'],mpileup,sample_list,cfg['reference']['FASTA'],variantcalling_log,workdir)
		varscan_vcf = tools.Concat_VarScan_vcf(varscan_snp_vcf,varscan_indel_vcf,variantcalling_log)
		varscan_vcf = tools.header_fix(dirs['script']+'merge_vcf.py',varscan_vcf,'V',variantcalling_log)
		varscan_vcf = tools.vcf_norm(cfg['tools']['BEDTOOLS']['path'],varscan_vcf,cfg['reference']['FASTA'],variantcalling_log)

		merge_vcf = tools.merge_vcf(dirs['script']+'merge_vcf.py',name,gatk_vcf,freeb_vcf,varscan_vcf,variantcalling_log,workdir)

		new_samplesheet.write('\t',join([gatk_vcf,freeb_vcf,varscan_vcf]))
		new_samplesheet.write(merge_vcf)
		preprocessing_log.close()
		new_samplesheet.close()

	if 'F' in workflow:
		pass

	if 'E' in workflow:
		pass



def Pipeline_Germline_Singlesample(workflow,samplesheet,design,panel,dirs,cfg,opts,target_list,target_bed,transcripts_list):

	status = subprocess.call("cat "+ dirs['logo']+'/logo_cmg.txt', shell=True)

	if 'A' in workflow:

		samples = f.ReadSampleSheet(samplesheet,'Germline',panel,'Alignment')
		alignment_log = open(dirs['log'] + '/Alignment.log','w+')
		samplesheet = dirs['log'] + '/preprocessing.samplesheet'
		new_samplesheet = open(samplesheet,'w+')

		f.makedirs([dirs['alignment']])
		workdir = dirs['alignment']

		for sample in samples.keys():
			if panel == 'CustomHPHS':
				sample_name,fq1,fq2,fqI2 = samples[sample][0]
				bam = Alignment(panel,sample_name,fq1,fq2,fqI2,dirs,cfg,opts,alignment_log,workdir)
			else:
				sample_name,fq1,fq2 = samples[sample][0]
				bam = Alignment(panel,sample_name,fq1,fq2,None,dirs,cfg,opts,alignment_log,workdir)

			new_samplesheet.write(sample_name+'\t'+bam+'\n')
		alignment_log.close()
		new_samplesheet.close()

	if 'R' in workflow or 'M' in workflow or 'I' in workflow or 'B' in workflow:

		samples = f.ReadSampleSheet(samplesheet,'Germline',panel,'Preprocessing')
		preprocessing_log = open(dirs['log'] + '/Preprocessing.log','w+')
		samplesheet = dirs['log'] + '/variantcalling.samplesheet'
		new_samplesheet = open(samplesheet,'w+')

		f.makedirs([dirs['preprocessing']])
		workdir = dirs['preprocessing']

		for sample in samples.keys():
			
			sample_name,bam = samples[sample][0]
			bam = Preprocessing(panel,workflow,target_list,sample_name,bam,dirs,cfg,opts,preprocessing_log,workdir)

			new_samplesheet.write(sample_name+'\t'+bam+'\n')

		preprocessing_log.close()
		new_samplesheet.close()
	
	if 'V' in workflow:

		pass

	if 'F' in workflow:
		pass

	if 'E' in workflow:
		pass

def Pipeline_Somatic_Case_Control(workflow,samplesheet,design,panel,dirs,cfg,opts,target_list,target_bed,transcripts_list):

	status = subprocess.call("cat "+ dirs['logo']+'/logo_cmg.txt', shell=True)

	if 'A' in workflow:

		samples = f.ReadSampleSheet(samplesheet,'Somatic_Case_Control',panel,'Alignment')
		alignment_log = open(dirs['log'] + '/Alignment.log','w+')
		samplesheet = dirs['log'] + '/preprocessing.samplesheet'
		new_samplesheet = open(samplesheet,'w+')

		f.makedirs([dirs['alignment']])
		workdir = dirs['alignment']

		for sample in samples.keys():
			if panel == 'CustomHPHS':
				sample_name,fq1,fq2,fqI2 = samples[sample][0]
				bam = Alignment(panel,sample_name,fq1,fq2,fqI2,dirs,cfg,opts,alignment_log,workdir)
			else:
				sample_name,fq1,fq2 = samples[sample][0]
				bam = Alignment(panel,sample_name,fq1,fq2,None,dirs,cfg,opts,alignment_log,workdir)

			new_samplesheet.write(sample_name+'\t'+bam+'\n')
		alignment_log.close()
		new_samplesheet.close()

	if 'R' in workflow or 'M' in workflow or 'I' in workflow or 'B' in workflow:

		samples = f.ReadSampleSheet(samplesheet,'Somatic_Case_Control',panel,'Preprocessing')
		preprocessing_log = open(dirs['log'] + '/Preprocessing.log','w+')
		samplesheet = dirs['log'] + '/variantcalling.samplesheet'
		new_samplesheet = open(samplesheet,'w+')

		f.makedirs([dirs['preprocessing']])
		workdir = dirs['preprocessing']

		for sample in samples.keys():
			
			sample_name,bam = samples[sample][0]
			bam = Preprocessing(panel,workflow,target_list,sample_name,bam,dirs,cfg,opts,preprocessing_log,workdir)

			new_samplesheet.write(sample_name+'\t'+bam+'\n')

		preprocessing_log.close()
		new_samplesheet.close()
	
	if 'V' in workflow:

		pass

	if 'F' in workflow:
		pass

	if 'E' in workflow:
		pass

def Pipeline_Somatic(workflow,samplesheet,design,panel,dirs,cfg,opts,target_list,target_bed,transcripts_list):

	status = subprocess.call("cat "+ dirs['logo']+'/logo_cmg.txt', shell=True)

	if 'A' in workflow:

		samples = f.ReadSampleSheet(samplesheet,'Somatic',panel,'Alignment')
		alignment_log = open(dirs['log'] + '/Alignment.log','w+')
		samplesheet = dirs['log'] + '/preprocessing.samplesheet'
		new_samplesheet = open(samplesheet,'w+')

		f.makedirs([dirs['alignment']])
		workdir = dirs['alignment']

		for sample in samples.keys():
			if panel == 'CustomHPHS':
				sample_name,fq1,fq2,fqI2 = samples[sample][0]
				bam = Alignment(panel,sample_name,fq1,fq2,fqI2,dirs,cfg,opts,alignment_log,workdir)
			else:
				sample_name,fq1,fq2 = samples[sample][0]
				bam = Alignment(panel,sample_name,fq1,fq2,None,dirs,cfg,opts,alignment_log,workdir)

			new_samplesheet.write(sample_name+'\t'+bam+'\n')
		alignment_log.close()
		new_samplesheet.close()

	if 'R' in workflow or 'M' in workflow or 'I' in workflow or 'B' in workflow:

		samples = f.ReadSampleSheet(samplesheet,'Somatic',panel,'Preprocessing')
		preprocessing_log = open(dirs['log'] + '/Preprocessing.log','w+')
		samplesheet = dirs['log'] + '/variantcalling.samplesheet'
		new_samplesheet = open(samplesheet,'w+')

		f.makedirs([dirs['preprocessing']])
		workdir = dirs['preprocessing']

		for sample in samples.keys():
			
			sample_name,bam = samples[sample][0]
			bam = Preprocessing(panel,workflow,target_list,sample_name,bam,dirs,cfg,opts,preprocessing_log,workdir)

			new_samplesheet.write(sample_name+'\t'+bam+'\n')

		preprocessing_log.close()
		new_samplesheet.close()
	
	if 'V' in workflow:

		pass

	if 'F' in workflow:
		pass

	if 'E' in workflow:
		pass