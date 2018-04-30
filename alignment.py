import subprocess
import Functions as f
import os

def ReadSampleSheet(samplesheet,analysis,panel):
	ssheet = dict()
	
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
					ssheet[sample_name] = [[gsample_name,gfq1,gfq2,gfqI2],[ssample_name,sfq1,sfq2,sfqI2]]
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
					ssheet[sample_name] = [[gsample_name,gfq1,gfq2],[ssample_name,sfq1,sfq2]]
			
	return ssheet

def Bwa_mem(path,threads,sample_name,fastq1,fastq2,ref,log,workdir):
	
	sam = workdir +'/'+sample_name+'.sam'
	sam_file=open(sam,'w+')

	print '- Alignment using BWA mem: '+sample_name
	
	if fastq2 == None:
		args = [path,'mem',opts.ref,fastq1,'-t',threads]
	else:
		args = [path,'mem',opts.ref,fastq1,fastq2,'-t',threads]
	
	success = subprocess.call(args, stdout=sam_file,stderr=log)
	if not success:
		return sam
	else:
		prRed('Error in Alignment. Check log file.')
		exit(1)

def SamFormatConverter(path,ram,sam,log,workdir):
	print '- From Sam to Bam'
	bam= workdir + '\t'+ '.'.join(sam.split('/')[-1].split('.')[:-1])+'.bam'
	args = ['java','-Xmx',ram,'-jar',path,'SamFormatConverter','I='+sam,'O='+bam]
	success = subprocess.call(args,stdout=log,stderr=log)
	if not success:
		return bam
	else:
		prRed('Error in Conversion. Check log file.')
		exit(1)

def SortSam(path,ram,bam,log,workdir):
	print '- Bam sorting'
	sort = workdir +'/'+ '.'.join(bam.split('/')[-1].split('.')[:-1])+'.sort.bam'
	args = ['java','-Xmx',ram,'-jar',path,'SortSam','I='+bam,'O='+sort,'SORT_ORDER=coordinate']
	success = subprocess.call(args,stdout=log,stderr=log)
	if not success:
		return sort
	else:
		prRed('Error in Sorting. Check log file.')
		exit(1)

def Index_bam(path,ram,bam,log):
	print '- Bam indexing'
	bai = bam+'.bai'
	if os.path.exists(bai):
		status = subprocess.call("rm "+ bai , shell=True)
	args = ['java','-Xmx',ram,'-jar',path,'BuildBamIndex','I='+bam,'O='+bai,'VALIDATION_STRINGENCY=LENIENT']
	success = subprocess.call(args,stdout=log,stderr=log)
	if not success:
		pass
	else:
		prRed('Error in Indexing. Check log file.')
		exit(1)

def SureCallTrimmer(path,ram,fq1,fq2,tag,outdir,log,workdir):
	print '- Adapters trimming'
	args = ['java','-Xmx',ram,'-jar',path,'-fq1',fq1,'-fq2',fq2,'-'+tag,'-out_loc',outdir]
	success = subprocess.call(args,stdout=log,stderr=log)

	fq1_name= fq1.split('/')[-1]
	fq2_name= fq2.split('/')[-1]

	for file in os.listdir(outdir):
		if fq1_name in file:
			trimmed_fq1 = outdir+'/'+file
		elif fq2_name in file:
			trimmed_fq2 = outdir+'/'+file

	if not success:
		return trimmed_fq1,trimmed_fq2
	else:
		prRed('Error in Trimming of adapters. Check log file.')
		exit(1)
	
def LocatIt(path,ram,bam,fqI2,log,workdir):
	print "- Merging haloplex molecular barcodes"
	outbam = workdir +'/'+ '.'.join(sam.split('/')[-1].split('.')[:-1] + ['MCBmerged'] +['bam'])
	args = ['java','-Xmx',ram,'-jar',path,'-U','-IB','-OB','-i','-o',outbam,bam,fqI2]
	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		return outbam
	else:
		prRed('Error in merging molucular barcodes. Check log file.')
		exit(1)

def Alignment(panel,sample_name,fq1,fq2=None,fqI2=None,dirs,cfg,opts,log):

	print 'ALIGNMENT sample: '+sample_name
	f.makedirs([dirs['alignment']])

	workdir = dirs['alignment']

	if panel == 'CustomSSQXT':

		fq1,fq2 = SureCallTrimmer(cfg['tools']['SURECALLTRIMMER']['path'],cfg['tools']['SURECALLTRIMMER']['ram'],fq1,fq2,'-hs',dirs['workdir']+'/TRIM_FASTQ',log)
		sam = Bwa_mem(cfg['tools']['BWA']['path'],cfg['tools']['BWA']['threads'],sample_name,fq1,fq2,cfg['reference']['FASTA'],log,workdir)
		bam = SamFormatConverter(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sam,log,workdir)
		mbbam = LocatIt(cfg['tools']['LOCATIT']['path'],cfg['tools']['LOCATIT']['ram'],bam,fqI2,log,workdir)
		sort = SortSam(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],mbbam,log,workdir)
		Index_bam(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sort,log)

		#f.Delete([sam,bam,mbbam])
		f.Move([sam,bam,mbbam],dirs['delete'])
	elif panel == 'CustomHPHS':

		fq1,fq2 = SureCallTrimmer(cfg['tools']['SURECALLTRIMMER']['path'],cfg['tools']['SURECALLTRIMMER']['ram'],fq1,fq2,'-qxt',dirs['workdir']+'/TRIM_FASTQ',log,)
		sam = Bwa_mem(cfg['tools']['BWA']['path'],cfg['tools']['BWA']['threads'],sample_name,fq1,fq2,cfg['reference']['FASTA'],log,workdir)
		bam = SamFormatConverter(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sam,log,workdir)
		sort = SortSam(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],bam,log,workdir)
		Index_bam(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sort,log,workdir)

		#f.Delete([sam,bam])
		f.Move([sam,bam],dirs['delete'])
	else:

		sam = Bwa_mem(cfg['tools']['BWA']['path'],cfg['tools']['BWA']['threads'],sample_name,fq1,fq2,cfg['reference']['FASTA'],log,workdir)
		bam = SamFormatConverter(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sam,log,workdir)
		sort = SortSam(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],bam,log,workdir)
		Index_bam(cfg['tools']['PICARD']['path'],cfg['tools']['PICARD']['ram'],sort,log,workdir)

		#f.Delete([sam,bam])
		f.Move([sam,bam],dirs['delete'])
	
	print 'ALIGNMENT sample: '+sample_name+ ' Done.'
	return sort