import subprocess
import os
import tools
import datetime
import functions as f

def Alignment(panel,sample_name,fq1,fq2,fqI2,dirs,cfg,opts,log,workdir):
	start_time = datetime.datetime.now()
	#print('Sample: '+sample_name
	f.makedirs([dirs['alignment']])

	workdir = dirs['alignment']

	if panel == 'CustomHPHS':
		f.makedirs([workdir+'/TRIM_FASTQ'])
		tr_fq1,tr_fq2 = tools.SureCallTrimmer(cfg['tools']['SURECALLTRIMMER']['path'],cfg['tools']['SURECALLTRIMMER']['ram'],fq1,fq2,'hs',workdir+'/TRIM_FASTQ',log,workdir)
		sam = tools.Bwa_mem(cfg['tools']['BWA']['path'],cfg['tools']['BWA']['threads'],sample_name,tr_fq1,tr_fq2,cfg['reference']['FASTA'],log,workdir)
		bam = tools.SamFormatConverter(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sam,log,workdir)
		mbbam = tools.LocatIt(cfg['tools']['LOCATIT']['path'],cfg['tools']['LOCATIT']['ram'],bam,fqI2,log,workdir)
		sort = tools.SortSam(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],mbbam,log,workdir)
		tools.BuildBamIndex(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sort,log)

		#f.Delete([sam,bam,mbbam])
		f.Move([sam,bam,mbbam],dirs['delete'])
		f.Copy([fq1,fq2,fqI2],dirs['storage'])
	elif panel == 'CustomSSQXT':
		f.makedirs([workdir+'/TRIM_FASTQ'])
		tr_fq1,tr_fq2 = tools.SureCallTrimmer(cfg['tools']['SURECALLTRIMMER']['path'],cfg['tools']['SURECALLTRIMMER']['ram'],fq1,fq2,'qxt',workdir+'/TRIM_FASTQ',log,workdir)
		sam = tools.Bwa_mem(cfg['tools']['BWA']['path'],cfg['tools']['BWA']['threads'],sample_name,tr_fq1,tr_fq2,cfg['reference']['FASTA'],log,workdir)
		bam = tools.SamFormatConverter(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sam,log,workdir)
		sort = tools.SortSam(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],bam,log,workdir)
		tools.BuildBamIndex(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sort,log)

		#f.Delete([sam,bam])
		f.Move([sam,bam],dirs['delete'])
		f.Copy([fq1,fq2],dirs['storage'])
	else:
		sam = tools.Bwa_mem(cfg['tools']['BWA']['path'],cfg['tools']['BWA']['threads'],sample_name,fq1,fq2,cfg['reference']['FASTA'],log,workdir)
		bam = tools.SamFormatConverter(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sam,log,workdir)
		sort = tools.SortSam(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],bam,log,workdir)
		tools.BuildBamIndex(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sort,log)

		#f.Delete([sam,bam])
		f.Move([sam,bam],dirs['delete'])
		f.Copy([fq1,fq2],dirs['storage'])
	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)

	print('Sample: '+sample_name+ ' -> Done in: %d min, %d sec' % elapsed_time)
	return sort

def Preprocessing(panel,workflow,target_list,sample_name,bam,dirs,cfg,opts,log,workdir):

	#print('Sample: ' + sample_name)
	start_time = datetime.datetime.now()

	if 'R' in workflow:

		ARRG_bam = tools.AddOrReplaceReadGroups(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],bam,sample_name,panel,opts.run_id,log,workdir)
		tools.BuildBamIndex(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],ARRG_bam,log)
		f.Move([bam,bam+'.bai'],dirs['delete'])
		bam = ARRG_bam

	if 'M' in workflow:

		MD_bam = tools.MarkDuplicates(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],bam,log,workdir)
		tools.BuildBamIndex(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],MD_bam,log)
		#f.Delete([sam,bam])
		f.Move([bam,bam+'.bai'],dirs['delete'])
		bam = MD_bam

	if 'I' in workflow:

		IR_bam = tools.IndelRealigner(cfg['tools']['GATK']['path'],cfg['tools']['GATK']['ram'],bam,cfg['reference']['FASTA'],cfg['database']['MILLS'],target_list,log,workdir)
		tools.BuildBamIndex(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],IR_bam,log)
		#f.Delete([sam,bam])
		f.Move([bam,bam+'.bai'],dirs['delete'])
		bam = IR_bam

	if 'B' in workflow:

		BR_bam = tools.BaseRecalibrator(cfg['tools']['GATK']['path'],cfg['tools']['GATK']['ram'],bam,cfg['reference']['FASTA'],cfg['database']['DBSNP'],cfg['database']['MILLS'],target_list,log,workdir)
		#f.Delete([sam,bam])
		f.Move([bam,bam+'.bai'],dirs['delete'])
		bam = BR_bam

	tools.BuildBamIndex(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],bam,log)

	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)
	print('Sample: '+sample_name+ ' -> Done in: %d min, %d sec' % elapsed_time)
	return bam
	#f.Delete([sam,bam])

def Pipeline_Germline_Multisample(workflow,samplesheet,design,panel,dirs,cfg,opts,target_list,target_bed,transcripts_list,cnv_target_list,cnv_target_bed,cnv_ref_ploidy,cnv_ref_calls):

	#status = subprocess.call("cat "+ dirs['logo']+'/logo_cmg.txt', shell=True)
	print('\n\nPOWERCALL v3.2.0 IS STARTING')
	print('"The coolest Pipeline I have ever seen" (Cit. Anonymous Fan)')
	print('Powered by Mario Urtis and CMGCV, IRCCS Fondazione Policlinico S.Matteo, Pavia, (IT)\n\n')

	if 'A' in workflow:

		print("ALIGNMENT:")
		if panel == 'CustomSSQXT' or panel == 'CustomHPHS': print("-Adapters trimming")
		print("-BWA mem")
		print("-SamFormatConverter")
		if panel == 'CustomHPHS': print("-Merge haloplex molecular barcodes")
		print("-Sort Bam\n")

		samples = f.ReadSampleSheet(samplesheet,'Germline',panel,'Alignment')
		alignment_log = open(dirs['log'] + '/Alignment.log','w+')
		samplesheet = dirs['log'] + '/Preprocessing.samplesheet'
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
		print("\nPREPROCESSING:")

		if "R" in workflow: print("-AddOrReplaceReadGroups")
		if "M" in workflow: print("-MarkDuplicates")
		if "I" in workflow: print("-IndelRealigner")
		if "B" in workflow: print("-BaseRecalibrator")
		print("")
		samples = f.ReadSampleSheet(samplesheet,'Germline',panel,'Preprocessing')
		preprocessing_log = open(dirs['log'] + '/Preprocessing.log','w+')
		samplesheet = dirs['log'] + '/Variantcalling.samplesheet'
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
		print("\nVARIANT CALLING")

		samples = f.ReadSampleSheet(samplesheet,'Germline',panel,'VariantCalling')
		variantcalling_log = open(dirs['log'] + '/VariantCalling.log','w+')
		samplesheet = dirs['log'] + '/Featuresextraction.samplesheet'
		CNVsamplesheet = open(dirs['log'] + '/CopyNumber.samplesheet','w+')
		
		new_samplesheet = open(samplesheet,'w+')

		workdir = dirs['variantcalling']
		f.makedirs([dirs['variantcalling'],dirs['gvcf']])

		bam_list = workdir+'/bams.list'
		sample_list = workdir+'/samples.list'
		gvcf_list = workdir+'/gvcf.list'

		bam_l = open(bam_list,'w+')
		sample_l = open(sample_list,'w+')
		gvcf_l = open(gvcf_list,'w+')
		bam_array = []
		gvcf_array = []
		to_ss = []

		for sample in samples.keys():
			sample_name,bam = samples[sample][0]
			gvcf = tools.HaplotypeCaller(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],bam,sample_name,cfg['reference']['FASTA'],target_list,variantcalling_log,dirs['gvcf'])
			bam_array+=[[bam,sample_name]]
			gvcf_array += [gvcf]
			bam_l.write(bam+'\n')
			sample_l.write(sample_name+'\n')
			gvcf_l.write(gvcf+'\n')
		
		f.Copy(gvcf_array,dirs['storage'])

		bam_l.close()
		sample_l.close()
		gvcf_l.close()

		name = opts.run_id

		gvcf_cohort = tools.CombineGVCFs(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],gvcf_array,name,cfg['reference']['FASTA'],target_list,variantcalling_log,dirs['gvcf'])
		gatk_vcf = tools.GenotypeGVCFs(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],gvcf_cohort,name,cfg['reference']['FASTA'],target_bed,variantcalling_log,workdir)
		gatk_vcf = tools.header_fix(dirs['script']+'header_fix.py',gatk_vcf,'G',variantcalling_log)
		gatk_vcf = tools.vcf_norm(cfg['tools']['BCFTOOLS']['path'],gatk_vcf,cfg['reference']['FASTA'],variantcalling_log)
		gatk_vcf = tools.HardFilter(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],gatk_vcf,name,cfg['reference']['FASTA'],target_bed,variantcalling_log,workdir)
		to_ss += [gatk_vcf]

		if not opts.GATKBestPractices:
		
			freeb_vcf = tools.FreeBayes(cfg['variantcaller']['FREEBAYES']['path'],name,None,bam_list,cfg['reference']['FASTA'],target_bed,variantcalling_log,workdir)
			freeb_vcf = tools.header_fix(dirs['script']+'header_fix.py',freeb_vcf,'F',variantcalling_log)
			freeb_vcf = tools.vcf_norm(cfg['tools']['BCFTOOLS']['path'],freeb_vcf,cfg['reference']['FASTA'],variantcalling_log)
			freeb_vcf = tools.HardFilter(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],freeb_vcf,name,cfg['reference']['FASTA'],target_bed,variantcalling_log,workdir)
			to_ss += [freeb_vcf]
			
			# mpileup = tools.mpileup(name,cfg['reference']['FASTA'],None,None,bam_list,target_bed,variantcalling_log,workdir)
			# varscan_snp_vcf = tools.VarScan_mpileup2snp(cfg['variantcaller']['VARSCAN']['path'],cfg['variantcaller']['VARSCAN']['ram'],mpileup,sample_list,cfg['reference']['FASTA'],target_bed,variantcalling_log,workdir)
			# varscan_indel_vcf = tools.VarScan_mpileup2indel(cfg['variantcaller']['VARSCAN']['path'],cfg['variantcaller']['VARSCAN']['ram'],mpileup,sample_list,cfg['reference']['FASTA'],target_bed,variantcalling_log,workdir)
			# varscan_vcf = tools.Concat_VarScan_vcf(varscan_snp_vcf,varscan_indel_vcf,variantcalling_log)
			# varscan_vcf = tools.header_fix(dirs['script']+'header_fix.py',varscan_vcf,'V',variantcalling_log)
			# varscan_vcf = tools.vcf_norm(cfg['tools']['BCFTOOLS']['path'],varscan_vcf,cfg['reference']['FASTA'],variantcalling_log)
			# to_ss += [varscan_vcf]

			merge_vcf = tools.merge_vcf(dirs['script']+'merge_vcf.py',name,gatk_vcf,freeb_vcf,None,variantcalling_log,workdir)
			to_ss = [merge_vcf]

		for bam,sample_name in bam_array:
			bai= bam+'.bai'
			storage_bam = dirs['storage'] + '/'+sample_name+'.bam'
			storage_bai = dirs['storage'] + '/'+sample_name+'.bam.bai'
			status = subprocess.call("cp " + bam + ' '+storage_bam , shell=True)
			status = subprocess.call("cp " + bai + ' '+storage_bai , shell=True)
			CNVsamplesheet.write('\t'.join([sample_name,storage_bam])+'\n')

		f.Copy(to_ss,dirs['storage'])
		new_samplesheet.write('\t'.join([name]+to_ss))
		variantcalling_log.close()
		CNVsamplesheet.close()
		new_samplesheet.close()

	if 'F' in workflow:

		print("\nFEATURES EXTRACTION")

		samples = f.ReadSampleSheet(samplesheet,'Germline',panel,'Featuresextraction')
		featuresextraction_log = open(dirs['log'] + '/FeaturesExtraction.log','w+')
		samplesheet = dirs['log'] + '/Annotation.samplesheet'
		new_samplesheet = open(samplesheet,'w+')
		f.makedirs([dirs['featsextract'],dirs['variantcalling']])
		workdir = dirs['variantcalling']

		if opts.GATKBestPractices:
			for sample in samples.keys():
				name, gatk_vcf = samples[sample][0]
				workdir = dirs['featsextract']
				toannotate_vcf, samples_list = tools.features_extractor(dirs['script']+'features_extractor_GATK.py',workdir,gatk_vcf,None,None,None,cfg['files']['LISTAFEATURES_GERMLINE'],dirs['gvcf'],design,featuresextraction_log,workdir)
		else:
			for sample in samples.keys():
				name, merge_vcf = samples[sample][0]
				bam_list = workdir+'/bams.list'
				#ieva_vcf = tools.iEVa(cfg['tools']['iEVA']['path'],merge_vcf,cfg['reference']['FASTA'],bam_list,log,workdir)
				workdir = dirs['featsextract']
				toannotate_vcf, samples_list = tools.features_extractor(dirs['script']+'features_extractor.py',workdir,None,None,None,merge_vcf,cfg['files']['LISTAFEATURES_GERMLINE'],dirs['gvcf'],design,featuresextraction_log,workdir)
			
		new_samplesheet.write('\t'.join([name,toannotate_vcf,samples_list]))
		new_samplesheet.close()
		featuresextraction_log.close()

	if 'E' in workflow:

		print("\nANNOTATION")

		samples = f.ReadSampleSheet(samplesheet,'Germline',panel,'Annotation')
		annotation_log = open(dirs['log'] + '/Annotation.log','w+')
		workdir = dirs['annotation']
		f.makedirs([dirs['annotation']])
		sample_list = workdir+'/samples.list'
		for sample in samples.keys():
			name,vcf,tsvfile= samples[sample][0]
			annotated_vcf = tools.Vep(cfg['annotation']['VEP']['path'],cfg['annotation']['VEP']['fork'],name,cfg['reference']['FASTA'],transcripts_list,cfg['annotation']['VEP']['af'],cfg['annotation']['VEP']['plugins'],vcf,annotation_log,workdir)
			workdir = dirs['out']
			for samp in open(tsvfile,'r'):
				samp = samp.rstrip()
				name = samp.split('/')[-1].split('.')[0]
				out_tsv = tools.add_Annotation(dirs['script']+'annotation_extractor.py',name,annotated_vcf,samp,cfg['files']['ANN_LIST_GERMLINE'],transcripts_list,annotation_log,workdir)
				f.Copy([out_tsv],dirs['storage'])
		
		samplesheet = dirs['log'] + '/CopyNumber.samplesheet'

	if 'C' in workflow:

		print("\nCNV CALLING")		
		#samplesheet = dirs['log'] + '/CopyNumber.samplesheet'
		samples = f.ReadSampleSheet(samplesheet,'Germline',panel,'CNV')
		cnv_log = open(dirs['log'] + '/CopyNumber.log','w+')
		name = opts.run_id
		f.makedirs([dirs['CNV'], dirs['CNV_HDF5'], dirs['CNV_PLOIDY'], dirs['CNV_CALLS']])
		workdir = dirs['CNV']
		bam_list = workdir+'/bams.list'
		bam_l = open(bam_list,'w+')
		sample_index = 0
		hdf5_array = []
		cnv_vcf_array = []
		for sample in samples.keys():
			sample_name,bam = samples[sample][0]
			hdf5 = tools.GATK_CollectReadCounts(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],bam,sample_name,cnv_target_list,cnv_log,dirs['CNV_HDF5'])
			hdf5_array += [[hdf5,sample_name,sample_index]]
			f.Copy([hdf5],dirs['storage'])
			bam_l.write(bam+'\n')
			sample_index += 1
		
		bam_l.close()

		if cnv_ref_ploidy == "":
			cnv_ref_ploidy = dirs['files'] +'/Prior.ploidy'
			print(cnv_ref_ploidy)

		GATK_sample_ploidy = tools.GATK_DetermineGermlineContigPloidy(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],hdf5_array,name,cnv_ref_ploidy,cnv_log,dirs['CNV_PLOIDY'])
		GATK_sample_calls, GATK_sample_model = tools.GATK_GermlineCNVCaller(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],hdf5_array,name,GATK_sample_ploidy,cnv_ref_calls,cnv_log,dirs['CNV_CALLS'])

		if cnv_ref_calls == "":
			cnv_ref_calls = GATK_sample_model

		for hdf5,sample_name,sample_index in hdf5_array:
			cnv_vcf, segments_vcf = tools.GATK_PostprocessGermlineCNVCalls(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],sample_name,sample_index,GATK_sample_calls,GATK_sample_ploidy,cnv_ref_calls,cnv_log,dirs['CNV'])
			cnv_vcf_array += [cnv_vcf]
			f.Copy([cnv_vcf, segments_vcf],dirs['storage'])
			f.Copy([cnv_vcf, segments_vcf],dirs['out'])
		

		GATK_CNV_report = tools.Make_CNV_report(dirs['files'] +'/'+panel +'_exons_info.txt', cnv_vcf_array, cnv_log, dirs['CNV'])
		f.Copy([GATK_CNV_report],dirs['out'])
		f.Copy([GATK_CNV_report],dirs['storage'])
			
		# print("-DECoN")
		# with f.working_directory(cfg['variantcaller']['DECON']['path']):
		# 	rdata = tools.Decon_ReadInBams(cfg['variantcaller']['DECON']['path'],bam_list,name,cnv_target_bed,cfg['reference']['FASTA'],cnv_log,dirs['CNV'])
		# 	tools.Decon_IdentifyFailures(cfg['variantcaller']['GATK']['path'],rdata,name,cfg['reference']['FASTA'],cnv_log,dirs['CNV'])
		# 	decon_calls = tools.Decon_makeCNVcalls(cfg['variantcaller']['GATK']['path'],rdata,name,cfg['reference']['FASTA'],cnv_log,dirs['CNV'])
		# 	f.Copy([decon_calls],dirs['storage'])
		# 	f.Copy([decon_calls],dirs['out'])

def Pipeline_Germline_Singlesample(workflow,samplesheet,design,panel,dirs,cfg,opts,target_list,target_bed,transcripts_list,cnv_target_list,cnv_target_bed,cnv_ref_ploidy,cnv_ref_calls):

	status = subprocess.call("cat "+ dirs['logo']+'/logo_cmg.txt', shell=True)

	if 'A' in workflow:

		samples = f.ReadSampleSheet(samplesheet,'Germline',panel,'Alignment')
		alignment_log = open(dirs['log'] + '/Alignment.log','w+')
		samplesheet = dirs['log'] + '/Preprocessing.samplesheet'
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
		samplesheet = dirs['log'] + '/Variantcalling.samplesheet'
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
		print("\nVARIANT CALLING")

		samples = f.ReadSampleSheet(samplesheet,'Germline',panel,'VariantCalling')
		variantcalling_log = open(dirs['log'] + '/VariantCalling.log','w+')
		samplesheet = dirs['log'] + '/Featuresextraction.samplesheet'
		
		new_samplesheet = open(samplesheet,'w+')

		workdir = dirs['variantcalling']
		f.makedirs([dirs['variantcalling'],dirs['gvcf']])

		bam_list = workdir+'/bams.list'
		sample_list = workdir+'/samples.list'
		gvcf_list = workdir+'/gvcf.list'

		bam_l = open(bam_list,'w+')
		sample_l = open(sample_list,'w+')
		gvcf_l = open(gvcf_list,'w+')
		bam_array = []
		gvcf_array = []
		to_ss = []

		for sample in samples.keys():
			sample_name,bam = samples[sample][0]
			gvcf = tools.HaplotypeCaller(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],bam,sample_name,cfg['reference']['FASTA'],target_list,variantcalling_log,dirs['gvcf'])
			gatk_vcf = tools.GenotypeGVCFs(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],gvcf,sample_name,cfg['reference']['FASTA'],target_bed,variantcalling_log,workdir)
			gatk_vcf = tools.header_fix(dirs['script']+'header_fix.py',gatk_vcf,'G',variantcalling_log)
			gatk_vcf = tools.vcf_norm(cfg['tools']['BCFTOOLS']['path'],gatk_vcf,cfg['reference']['FASTA'],variantcalling_log)

			bam_array+=[[bam,sample_name]]
			gvcf_array += [gvcf]
			bam_l.write(bam+'\n')
			sample_l.write(sample_name+'\n')
			gvcf_l.write(gvcf+'\n')
			f.Copy([gvcf],dirs['storage'])

		bam_l.close()
		sample_l.close()
		gvcf_l.close()

		name = opts.run_id


	if 'F' in workflow:
		pass

	if 'E' in workflow:
		pass

def Pipeline_Somatic_Case_Control(workflow,samplesheet,design,panel,dirs,cfg,opts,target_list,target_bed,transcripts_list,cnv_target_list,cnv_target_bed,cnv_ref_ploidy,cnv_ref_calls):

	status = subprocess.call("cat "+ dirs['logo']+'/logo_cmg.txt', shell=True)

	if 'A' in workflow:

		print("ALIGNMENT:")
		if panel == 'CustomSSQXT' or panel == 'CustomHPHS': print("-Adapters trimming")
		print("-BWA mem")
		print("-SamFormatConverter")
		if panel == 'CustomHPHS': print("-Merge haloplex molecular barcodes")
		print("-Sort Bam\n")

		samples = f.ReadSampleSheet(samplesheet,'Somatic_Case_Control',panel,'Alignment')
		alignment_log = open(dirs['log'] + '/Alignment.log','w+')
		samplesheet = dirs['log'] + '/Preprocessing.samplesheet'
		new_samplesheet = open(samplesheet,'w+')

		f.makedirs([dirs['alignment']])
		workdir = dirs['alignment']

		for sample in samples.keys():
			if panel == 'CustomHPHS':
				gsample_name,gfq1,gfq2,gfqI2 = samples[sample][0]
				gbam = Alignment(panel,gsample_name,gfq1,gfq2,gfqI2,dirs,cfg,opts,alignment_log,workdir)

				ssample_name,sfq1,sfq2,sfqI2 = samples[sample][1]
				sbam = Alignment(panel,ssample_name,sfq1,sfq2,sfqI2,dirs,cfg,opts,alignment_log,workdir)

			else:
				gsample_name,gfq1,gfq2 = samples[sample][0]
				gbam = Alignment(panel,gsample_name,gfq1,gfq2,None,dirs,cfg,opts,alignment_log,workdir)

				ssample_name,sfq1,sfq2 = samples[sample][1]
				sbam = Alignment(panel,ssample_name,sfq1,sfq2,None,dirs,cfg,opts,alignment_log,workdir)

			new_samplesheet.write(gsample_name+'\t'+gbam+'\t'+ssample_name+'\t'+sbam+'\n')
		alignment_log.close()
		new_samplesheet.close()

	if 'R' in workflow or 'M' in workflow or 'I' in workflow or 'B' in workflow:

		print("\nPREPROCESSING:")

		if "R" in workflow: print("-AddOrReplaceReadGroups")
		if "M" in workflow: print("-MarkDuplicates")
		if "I" in workflow: print("-IndelRealigner")
		if "B" in workflow: print("-BaseRecalibrator")

		samples = f.ReadSampleSheet(samplesheet,'Somatic_Case_Control',panel,'Preprocessing')
		preprocessing_log = open(dirs['log'] + '/Preprocessing.log','w+')
		samplesheet = dirs['log'] + '/Variantcalling.samplesheet'
		new_samplesheet = open(samplesheet,'w+')

		f.makedirs([dirs['preprocessing']])
		workdir = dirs['preprocessing']

		for sample in samples.keys():
			
			gsample_name,gbam = samples[sample][0]
			gbam = Preprocessing(panel,workflow,target_list,gsample_name,gbam,dirs,cfg,opts,preprocessing_log,workdir)

			ssample_name,sbam = samples[sample][1]
			sbam = Preprocessing(panel,workflow,target_list,ssample_name,sbam,dirs,cfg,opts,preprocessing_log,workdir)

			new_samplesheet.write(gsample_name+'\t'+gbam+'\t'+ssample_name+'\t'+sbam+'\n')

		preprocessing_log.close()
		new_samplesheet.close()
	
	if 'V' in workflow:
		print("\nVARIANT CALLING")

		samples = f.ReadSampleSheet(samplesheet,'Somatic_Case_Control',panel,'VariantCalling')
		variantcalling_log = open(dirs['log'] + '/VariantCalling.log','w+')
		samplesheet = dirs['log'] + '/Featuresextraction.samplesheet'
		
		new_samplesheet = open(samplesheet,'w+')

		workdir = dirs['variantcalling']
		f.makedirs([dirs['variantcalling']])


		for sample in samples.keys():

			gsample_name,gbam = samples[sample][0]
			ssample_name,sbam = samples[sample][1]
			print("-"+ssample_name)
			mutect_vcf = tools.Mutect2(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],gbam,sbam,gsample_name,ssample_name,cfg['reference']['FASTA'],target_bed,variantcalling_log,workdir)
			mutect_vcf = tools.vcf_norm(cfg['tools']['BCFTOOLS']['path'],mutect_vcf,cfg['reference']['FASTA'],variantcalling_log)

			vardict_vcf = tools.Vardict(cfg['variantcaller']['VARDICT']['path'],cfg['variantcaller']['VARDICT']['script_dir'],cfg['variantcaller']['VARDICT']['threads'],gbam,sbam,gsample_name,ssample_name,cfg['reference']['FASTA'],target_bed,variantcalling_log,workdir)
			vardict_vcf = tools.vcf_norm(cfg['tools']['BCFTOOLS']['path'],vardict_vcf,cfg['reference']['FASTA'],variantcalling_log)
			
			mpileup = tools.mpileup(ssample_name,cfg['reference']['FASTA'],gbam,sbam,None,target_bed,variantcalling_log,workdir)
			varscan_snp_vcf,varscan_indel_vcf = tools.VarScan_somatic(cfg['variantcaller']['VARSCAN']['path'],cfg['variantcaller']['VARSCAN']['ram'],None,mpileup,cfg['reference']['FASTA'],target_bed,variantcalling_log,workdir)
			varscan_vcf = tools.Concat_VarScan_vcf(varscan_snp_vcf,varscan_indel_vcf,variantcalling_log)
			varscan_vcf = tools.header_fix(dirs['script']+'header_fix.py',varscan_vcf,'V',variantcalling_log)
			varscan_vcf = tools.vcf_norm(cfg['tools']['BCFTOOLS']['path'],varscan_vcf,cfg['reference']['FASTA'],variantcalling_log)
			
			new_samplesheet.write('\t'.join([gsample_name,ssample_name,mutect_vcf,vardict_vcf,varscan_vcf]))

		variantcalling_log.close()
		new_samplesheet.close()

	if 'F' in workflow:
		print("\nFEATURES EXTRACTION")

		samples = f.ReadSampleSheet(samplesheet,'Somatic_Case_Control',panel,'Featuresextraction')
		featuresextraction_log = open(dirs['log'] + '/FeaturesExtraction.log','w+')
		samplesheet = dirs['log'] + '/Annotation.samplesheet'
		new_samplesheet = open(samplesheet,'w+')
		
		f.makedirs([dirs['featsextract']])
		workdir = dirs['variantcalling']

		for sample in samples.keys():
			gsample_name,ssample_name,mutect_vcf,vardict_vcf,varscan_vcf = samples[sample][0]
			bam_list = workdir+'/bams.list'
			#ieva_vcf = tools.iEVa(cfg['tools']['iEVA']['path'],merge_vcf,cfg['reference']['FASTA'],bam_list,log,workdir)
			workdir = dirs['featsextract']
			name = ssample_name+'_SomaticCC'
			#toannotate_vcf = '/home/jarvis/Scrivania/TEST/pipeline/test1/VARIANTCALLING/20180510_prova.merge.vcf'
			tsvfile,vcffile = tools.features_extractor_somatic(dirs['script']+'features_extractor_somatic.py',workdir+'/'+name,mutect_vcf,vardict_vcf,varscan_vcf,gsample_name,ssample_name,cfg['files']['LISTAFEATURES_SOMATIC'],target_bed,featuresextraction_log,workdir)
			new_samplesheet.write('\t'.join([name,vcffile,tsvfile]))
		new_samplesheet.close()
		featuresextraction_log.close()
	
	if 'E' in workflow:
		print("\nANNOTATION")

		samples = f.ReadSampleSheet(samplesheet,'Somatic_Case_Control',panel,'Annotation')
		annotation_log = open(dirs['log'] + '/Annotation.log','w+')
		#samplesheet = dirs['log'] + '/.samplesheet'
		#new_samplesheet = open(samplesheet,'w+')
		workdir = dirs['annotation']
		f.makedirs([dirs['annotation']])
		for sample in samples.keys():
			 #print(sample)
			name,vcf,tsvfile= samples[sample][0]
			annotated_vcf = tools.Vep(cfg['annotation']['VEP']['path'],cfg['annotation']['VEP']['fork'],name,cfg['reference']['FASTA'],transcripts_list,cfg['annotation']['VEP']['af'],cfg['annotation']['VEP']['plugins'],vcf,annotation_log,workdir)
			workdir = dirs['out']
			#annotated_vcf = '/home/jarvis/Scrivania/TEST/pipeline/test1/ANNOTATION/20180510_prova.VEP.vcf'
			tools.add_Annotation(dirs['script']+'annotation_extractor.py',name,annotated_vcf,tsvfile,cfg['files']['ANN_LIST_SOMATIC'],transcripts_list,annotation_log,workdir)
		pass

def Pipeline_Somatic(workflow,samplesheet,design,panel,dirs,cfg,opts,target_list,target_bed,transcripts_list,cnv_target_list,cnv_target_bed,cnv_ref_ploidy,cnv_ref_calls):

	#status = subprocess.call("cat "+ dirs['logo']+'/logo_cmg.txt', shell=True)
	print('\n\nPOWERCALL v3.2.0 IS STARTING')
	print('"The coolest Pipeline I have ever seen" (Cit. Anonymous Fan)')
	print('Powered by Mario Urtis and CMGCV, IRCCS Fondazione Policlinico S.Matteo, Pavia, (IT)\n\n')

	if 'A' in workflow:

		print("ALIGNMENT:")
		if panel == 'CustomSSQXT' or panel == 'CustomHPHS': print("-Adapters trimming")
		print("-BWA mem")
		print("-SamFormatConverter")
		if panel == 'CustomHPHS': print("-Merge haloplex molecular barcodes")
		print("-Sort Bam\n")

		samples = f.ReadSampleSheet(samplesheet,'Somatic',panel,'Alignment')
		alignment_log = open(dirs['log'] + '/Alignment.log','w+')
		samplesheet = dirs['log'] + '/Preprocessing.samplesheet'
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
		print("\nPREPROCESSING:")

		if "R" in workflow: print("-AddOrReplaceReadGroups")
		if "M" in workflow: print("-MarkDuplicates")
		if "I" in workflow: print("-IndelRealigner")
		if "B" in workflow: print("-BaseRecalibrator")
		print("")
		samples = f.ReadSampleSheet(samplesheet,'Somatic',panel,'Preprocessing')
		preprocessing_log = open(dirs['log'] + '/Preprocessing.log','w+')
		samplesheet = dirs['log'] + '/Variantcalling.samplesheet'
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

		print("\nVARIANT CALLING")

		samples = f.ReadSampleSheet(samplesheet,'Somatic',panel,'VariantCalling')
		variantcalling_log = open(dirs['log'] + '/VariantCalling.log','w+')
		samplesheet = dirs['log'] + '/Featuresextraction.samplesheet'
		
		new_samplesheet = open(samplesheet,'w+')

		workdir = dirs['variantcalling']
		f.makedirs([dirs['variantcalling']])

		bam_list = workdir+'/bams.list'
		sample_list = workdir+'/samples.list'

		bam_l = open(bam_list,'w+')
		sample_l = open(sample_list,'w+')
		bam_array = []

		for sample in samples.keys():
			ssample_name,sbam = samples[sample][0]
			mutect_vcf = tools.Mutect2(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],None,sbam,None,ssample_name,cfg['reference']['FASTA'],target_list,variantcalling_log,workdir)
			mutect_vcf = tools.header_fix(dirs['script']+'header_fix.py',mutect_vcf,'M',variantcalling_log)
			mutect_vcf = tools.vcf_norm(cfg['tools']['BCFTOOLS']['path'],mutect_vcf,cfg['reference']['FASTA'],variantcalling_log)
			#mutect_vcf = tools.LeftAlignAndTrimVariants(cfg['variantcaller']['GATK']['path'],mutect_vcf,cfg['reference']['FASTA'],variantcalling_log)
			bam_array+=[[sbam,ssample_name]]
			new_samplesheet.write('\t'.join([ssample_name,mutect_vcf])+'\n')

		# gatk_vcf = tools.GenotypeGVCFs(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],name,gvcf_list,cfg['reference']['FASTA'],target_bed,variantcalling_log,workdir)
		# gatk_vcf = tools.header_fix(dirs['script']+'header_fix.py',gatk_vcf,'G',variantcalling_log)
		# gatk_vcf = tools.vcf_norm(cfg['tools']['BCFTOOLS']['path'],gatk_vcf,cfg['reference']['FASTA'],variantcalling_log)
		
		# freeb_vcf = tools.FreeBayes(cfg['variantcaller']['FREEBAYES']['path'],name,None,bam_list,cfg['reference']['FASTA'],target_bed,variantcalling_log,workdir)
		# freeb_vcf = tools.header_fix(dirs['script']+'header_fix.py',freeb_vcf,'F',variantcalling_log)
		# freeb_vcf = tools.vcf_norm(cfg['tools']['BCFTOOLS']['path'],freeb_vcf,cfg['reference']['FASTA'],variantcalling_log)
		
		# mpileup = tools.mpileup(name,cfg['reference']['FASTA'],None,None,bam_list,target_bed,variantcalling_log,workdir)
		# varscan_snp_vcf = tools.VarScan_mpileup2snp(cfg['variantcaller']['VARSCAN']['path'],cfg['variantcaller']['VARSCAN']['ram'],mpileup,sample_list,cfg['reference']['FASTA'],target_bed,variantcalling_log,workdir)
		# varscan_indel_vcf = tools.VarScan_mpileup2indel(cfg['variantcaller']['VARSCAN']['path'],cfg['variantcaller']['VARSCAN']['ram'],mpileup,sample_list,cfg['reference']['FASTA'],target_bed,variantcalling_log,workdir)
		# varscan_vcf = tools.Concat_VarScan_vcf(varscan_snp_vcf,varscan_indel_vcf,variantcalling_log)
		# varscan_vcf = tools.header_fix(dirs['script']+'header_fix.py',varscan_vcf,'V',variantcalling_log)
		# varscan_vcf = tools.vcf_norm(cfg['tools']['BCFTOOLS']['path'],varscan_vcf,cfg['reference']['FASTA'],variantcalling_log)

		for bam,sample_name in bam_array:
			bai= bam+'.bai'
			storage_bam = dirs['storage'] + '/'+sample_name+'.bam'
			storage_bai = dirs['storage'] + '/'+sample_name+'.bam.bai'
			status = subprocess.call("mv " + bam + ' '+storage_bam , shell=True)
			status = subprocess.call("mv " + bai + ' '+storage_bai , shell=True)

		# f.Copy([gatk_vcf,freeb_vcf,varscan_vcf],dirs['storage'])
		
		
		variantcalling_log.close()
		new_samplesheet.close()

	if 'F' in workflow:
		print("\nFEATURES EXTRACTION")

		samples = f.ReadSampleSheet(samplesheet,'Somatic',panel,'Featuresextraction')
		featuresextraction_log = open(dirs['log'] + '/FeaturesExtraction.log','w+')
		samplesheet = dirs['log'] + '/Annotation.samplesheet'
		new_samplesheet = open(samplesheet,'w+')
		
		f.makedirs([dirs['featsextract']])
		workdir = dirs['variantcalling']

		for sample in samples.keys():
			ssample_name,mutect_vcf = samples[sample][0]
			bam_list = workdir+'/bams.list'
			#ieva_vcf = tools.iEVa(cfg['tools']['iEVA']['path'],merge_vcf,cfg['reference']['FASTA'],bam_list,log,workdir)
			workdir = dirs['featsextract']
			#toannotate_vcf = '/home/jarvis/Scrivania/TEST/pipeline/test1/VARIANTCALLING/20180510_prova.merge.vcf'
			tsvfile,vcffile = tools.features_extractor_somatic_MUTECT2(dirs['script']+'features_extractor_somatic_MUTECT2.py',workdir+'/'+ssample_name,mutect_vcf,ssample_name,cfg['files']['LISTAFEATURES_MUTECT2'],target_bed,featuresextraction_log,workdir)
			new_samplesheet.write('\t'.join([ssample_name,vcffile,tsvfile])+'\n')
		new_samplesheet.close()
		featuresextraction_log.close()
	
	if 'E' in workflow:
		print("\nANNOTATION")

		samples = f.ReadSampleSheet(samplesheet,'Somatic',panel,'Annotation')
		annotation_log = open(dirs['log'] + '/Annotation.log','w+')
		#samplesheet = dirs['log'] + '/.samplesheet'
		#new_samplesheet = open(samplesheet,'w+')
		workdir = dirs['annotation']
		f.makedirs([dirs['annotation']])
		for sample in samples.keys():
			name,vcf,tsvfile= samples[sample][0]
			annotated_vcf = tools.Vep(cfg['annotation']['VEP']['path'],cfg['annotation']['VEP']['fork'],name,cfg['reference']['FASTA'],transcripts_list,cfg['annotation']['VEP']['af'],cfg['annotation']['VEP']['plugins'],vcf,annotation_log,workdir)
			workdir = dirs['out']
			#annotated_vcf = '/home/jarvis/Scrivania/TEST/pipeline/test1/ANNOTATION/20180510_prova.VEP.vcf'
			tools.add_Annotation(dirs['script']+'annotation_extractor.py',name,annotated_vcf,tsvfile,cfg['files']['ANN_LIST_SOMATIC'],transcripts_list,annotation_log,workdir)
		pass