import subprocess
import Functions as f
import os


def ReadSampleSheet(samplesheet,analysis,panel):
	ssheet = dict()
	
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
				ssheet[sample_name] = [[gsample_name,gbam],[ssample_name,sbam]]

	return ssheet

def AddOrReplaceReadGroups(path,ram,bam,sample_name,panel,log,workdir):

	print "- Add/replace Read groups"
	outbam = workdir +'/'+ '.'.join(bam.split('/')[-1].split('.')[:-1] + ['ARRG'] +['bam'])
	args = ['java','-Xmx',ram,'-jar',path,'AddOrReplaceReadGroups','I='+bam,'O='+outbam,'RGID='+sample_name,'RGPL=ILLUMINA'+,'RGSM='sample_name+,'RGLB='+panel,'VALIDATION_STRINGENCY=LENIENT']
	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		return outbam
	else:
		prRed('Error in Adding or replacing Read groups. Check log file.')
		exit(1)

def BuildBamIndex(path,ram,bam,log,workdir):

	print '- Bam indexing'
	bai= bam+'.bai'
	if os.path.exists(bai):
		status = subprocess.call("rm "+ bai , shell=True)
	args = ['java','-Xmx',ram,'-jar',path,'BuildBamIndex','I='+bam,'O='+bai,'VALIDATION_STRINGENCY=LENIENT']
	success = subprocess.call(args,stdout=log,stderr=log)
	if not success:
		pass
	else:
		prRed('Error in Indexing. Check log file.')
		exit(1)


def MarkDuplicates(path,ram,bam,log,workdir):

	print "- Marking duplicates"
	outbam = workdir +'/'+ '.'.join(bam.split('/')[-1].split('.')[:-1] + ['Mark'] +['bam'])
	metrics_file = '.'.join(bam.split('/')[-1].split('.')[:-1] + ['MarkMetrics'] +['txt'])
	args = ['java','-Xmx',ram,'-jar',path,'MarkDuplicates','I='+bam,'O='+outbam,'METRICS_FILE='+metrics_file,'READ_NAME_REGEX=null','ASSUME_SORTED=true','VALIDATION_STRINGENCY=LENIENT']
	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		return outbam
	else:
		prRed('Error in Marking of duplicates. Check log file.')
		exit(1)

def IndelRealigner(path,ram,bam,reference,mills,target,log,workdir):

	print "- Indel realignment"
	outbam = workdir +'/'+ '.'.join(bam.split('/')[-1].split('.')[:-1] + ['Mark'] +['bam'])
	intervals = '.'.join(bam.split('/')[-1].split('.')[:-1] + ['IndelRealigner'] +['intervals'])
	args = ['java','-Xmx',ram,'-jar',path,'RealignerTargetCreator','-R',reference,'-I', bam,'.-o',intervals,'-known',mills,'-L',target]
	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		pass
	else:
		prRed('Error in target creation for indel realignment. Check log file.')
		exit(1)

	args = ['java','-Xmx',ram,'-jar',path,'IndelRealigner','-R',reference,'-I', bam,'-targetIntervals',intervals,'-known',mills,'-o',outbam]
	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		return outbam
	else:
		prRed('Error in indel realignment. Check log file.')
		exit(1)

def BaseRecalibrator(path,ram,bam,reference,dbsnp,mills,target,log,workdir):

	print "- Base quality score recalibration"
	outbam = workdir +'/'+ '.'.join(bam.split('/')[-1].split('.')[:-1] + ['BQSR'] +['bam'])
	table = workdir +'/'+ '.'.join(bam.split('/')[-1].split('.')[:-1] + ['BQSR'] +['table'])
	args = ['java','-Xmx',ram,'-jar',path,'BaseRecalibrator','-R',reference,'-I', bam,'.-o',table,'-knownSites',dbsnp,'-knownSites',mills,'-L',target]
	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		pass
	else:
		prRed('Error in table creation for Base quality score recalibration. Check log file.')
		exit(1)

	args = ['java','-Xmx',ram,'-jar',path,'PrintReads','-R',reference,'-I', bam,'.BQSR',table,'-L',target,'-o',outbam]
	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		return outbam
	else:
		prRed('Error in base quality score recalibration. Check log file.')
		exit(1)


def Preprocessing(panel,target,sample_name,bam,dirs,cfg,opts,log,target_list,target_bed,transcripts_list):

	print 'PREPROCESSING sample: ' + sample_name
	f.makedirs([dirs['preprocessing']])
	workdir = dirs['preprocessing']

	if 'R' in workflow:

		ARRG_bam = AddOrReplaceReadGroups(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],bam,sample_name,panel,log,workdir)
		BuildBamIndex((cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],ARRG_bam,log))
		bam = ARRG_bam

	if 'M' in workflow:

		MD_bam = MarkDuplicates(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],bam,log,workdir)
		#f.Delete([sam,bam])
		f.Move([bam],dirs['delete'])
		bam = MD_bam

	if 'I' in workflow:

		IR_bam = IndelRealigner(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],bam,cfg['reference']['FASTA'],cfg['database']['MILLS'],target_list,log,workdir)
		#f.Delete([sam,bam])
		f.Move([bam],dirs['delete'])
		bam = IR_bam

	if 'B' in workflow:

		BR_bam = BaseRecalibrator(cfg['variantcaller']['GATK']['path'],cfg['variantcaller']['GATK']['ram'],bam,cfg['reference']['FASTA'],cfg['database']['DBSNP'],cfg['database']['MILLS'],target_list,log,workdir)
		#f.Delete([sam,bam])
		f.Move([bam],dirs['delete'])
		bam = BR_bam

	BuildBamIndex((cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],bam,log))
	return bam
	#f.Delete([sam,bam])


