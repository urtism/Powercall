import subprocess
import Functions as f
import os
import datetime

###----------------------------------------------------------------------------------------------ALIGNMENT-------------------------------------------------------------------------------------


def Bwa_mem(path,threads,sample_name,fastq1,fastq2,reference,log,workdir):
	
	sam = workdir +'/'+sample_name+'.sam'
	sam_file=open(sam,'w+')

	for elem in [path,fastq1,fastq2,reference,workdir]:
		if os.path.isfile(elem) or os.path.isdir(elem):
			pass
		else: 
			print elem 
	#print '- BWA mem'
	
	if fastq2 == None:
		args = [path,'mem',reference,fastq1,'-t',threads]
	else:
		args = [path,'mem',reference,fastq1,fastq2,'-t',threads]
	
	success = subprocess.call(args, stdout=sam_file,stderr=log)
	if not success:
		return sam
	else:
		f.prRed('Error in Alignment. Check log file.')
		exit(1)

def SamFormatConverter(path,ram,sam,log,workdir):
	#print '- From Sam to Bam'
	bam= workdir + '/' + '.'.join(sam.split('/')[-1].split('.')[:-1])+'.bam'
	args = ['java','-Xmx'+ram,'-jar',path,'SamFormatConverter','I='+sam,'O='+bam]
	success = subprocess.call(args,stdout=log,stderr=log)
	if not success:
		return bam
	else:
		f.prRed('Error in Conversion. Check log file.')
		exit(1)

def SortSam(path,ram,bam,log,workdir):
	#print '- Bam sorting'
	sort = workdir +'/'+ '.'.join(bam.split('/')[-1].split('.')[:-1])+'.sort.bam'
	args = ['java','-Xmx'+ram,'-jar',path,'SortSam','I='+bam,'O='+sort,'SORT_ORDER=coordinate']
	success = subprocess.call(args,stdout=log,stderr=log)
	if not success:
		return sort
	else:
		f.prRed('Error in Sorting. Check log file.')
		exit(1)

def SureCallTrimmer(path,ram,fq1,fq2,tag,outdir,log,workdir):
	#print '- Adapters trimming'
	args = ['java','-Xmx'+ram,'-jar',path,'-fq1',fq1,'-fq2',fq2,'-'+tag,'-out_loc',outdir]
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
		f.prRed('Error in Trimming of adapters. Check log file.')
		exit(1)
	
def LocatIt(path,ram,bam,fqI2,log,workdir):
	#print "- Merging haloplex molecular barcodes"
	outbam = workdir +'/'+ '.'.join(sam.split('/')[-1].split('.')[:-1] + ['MCBmerged'] +['bam'])
	args = ['java','-Xmx'+ram,'-jar',path,'-U','-IB','-OB','-i','-o',outbam,bam,fqI2]
	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		return outbam
	else:
		f.prRed('Error in merging molucular barcodes. Check log file.')
		exit(1)

###----------------------------------------------------------------------------------------------PREPROCESSING-------------------------------------------------------------------------------------

def AddOrReplaceReadGroups(path,ram,bam,sample_name,panel,run,log,workdir):

	#print "- Add/replace Read groups"
	outbam = workdir +'/'+ '.'.join(bam.split('/')[-1].split('.')[:-1] + ['RG'] +['bam'])
	args = ['java','-Xmx'+ram,'-jar',path,'AddOrReplaceReadGroups','I='+bam,'O='+outbam,'RGID='+sample_name,'RGPL=ILLUMINA','RGSM='+sample_name,'RGLB='+panel,"RGPU="+run,'VALIDATION_STRINGENCY=LENIENT']
	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		return outbam
	else:
		f.prRed('Error in Adding or replacing Read groups. Check log file.')
		exit(1)

def BuildBamIndex(path,ram,bam,log):

	#print '- Bam indexing'
	bai= bam+'.bai'
	if os.path.exists(bai):
		status = subprocess.call("rm "+ bai , shell=True)
	args = ['java','-Xmx'+ram,'-jar',path,'BuildBamIndex','I='+bam,'O='+bai,'VALIDATION_STRINGENCY=LENIENT']
	success = subprocess.call(args,stdout=log,stderr=log)
	if not success:
		pass
	else:
		f.prRed('Error in Indexing. Check log file.')
		exit(1)


def MarkDuplicates(path,ram,bam,log,workdir):

	#print "- Marking duplicates"
	outbam = workdir +'/'+ '.'.join(bam.split('/')[-1].split('.')[:-1] + ['Mark'] +['bam'])
	metrics_file = workdir +'/'+ '.'.join(bam.split('/')[-1].split('.')[:-1] + ['MarkMetrics'] +['txt'])
	args = ['java','-Xmx'+ram,'-jar',path,'MarkDuplicates','I='+bam,'O='+outbam,'METRICS_FILE='+metrics_file,'READ_NAME_REGEX=null','ASSUME_SORTED=true','VALIDATION_STRINGENCY=LENIENT']
	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		return outbam
	else:
		f.prRed('Error in Marking of duplicates. Check log file.')
		exit(1)

def IndelRealigner(path,ram,bam,reference,mills,target,log,workdir):

	#print "- Indel realignment"
	outbam = workdir +'/'+ '.'.join(bam.split('/')[-1].split('.')[:-1] + ['IR'] +['bam'])
	intervals = workdir +'/'+ '.'.join(bam.split('/')[-1].split('.')[:-1] + ['IndelRealigner'] +['intervals'])
	args = ['java','-Xmx'+ram,'-jar',path,'-T','RealignerTargetCreator','-R',reference,'-I', bam,'-o',intervals,'-known',mills,'-L',target]
	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		pass
	else:
		f.prRed('Error in target creation for indel realignment. Check log file.')
		exit(1)

	args = ['java','-Xmx'+ram,'-jar',path,'-T','IndelRealigner','-R',reference,'-I', bam,'-targetIntervals',intervals,'-known',mills,'-o',outbam]
	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		return outbam
	else:
		f.prRed('Error in indel realignment. Check log file.')
		exit(1)

def BaseRecalibrator(path,ram,bam,reference,dbsnp,mills,target,log,workdir):

	#print "- Base quality score recalibration"
	outbam = workdir +'/'+ '.'.join(bam.split('/')[-1].split('.')[:-1] + ['BQSR'] +['bam'])
	table = workdir +'/'+ '.'.join(bam.split('/')[-1].split('.')[:-1] + ['BQSR'] +['table'])
	args = ['java','-Xmx'+ram,'-jar',path,'-T','BaseRecalibrator','-R',reference,'-I', bam,'-o',table,'-knownSites',dbsnp,'-knownSites',mills,'-L',target]
	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		pass
	else:
		f.prRed('Error in table creation for Base quality score recalibration. Check log file.')
		exit(1)

	args = ['java','-Xmx'+ram,'-jar',path,'-T','PrintReads','-R',reference,'-I', bam,'-BQSR',table,'-L',target,'-o',outbam]
	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		return outbam
	else:
		f.prRed('Error in base quality score recalibration. Check log file.')
		exit(1)



###----------------------------------------------------------------------------------------------VARIANTCALLING-------------------------------------------------------------------------------------

def mpileup(name,reference,gbam,sbam,bam_list,target,log,workdir):

	start_time = datetime.datetime.now()
	#print "- Mpileup"
	mpileup = workdir +'/'+name+'.mpileup'
	open_mpileup = open(mpileup,'w')
	if sbam != None:
		args = ['samtools','mpileup','-B','-q 1','-d','50000','-L','50000','-f',reference,gbam,sbam]
	elif gbam != None:
		args = ['samtools','mpileup','-B','-q 1','-d','50000','-L','50000','-f',reference,gbam]
	else:
		args = ['samtools','mpileup','-B','-q 1','-d','50000','-L','50000','-f',reference,'-b',bam_list]
	
	if target != None:
		args += ['-l',target]

	success = subprocess.call(args,stdout=open_mpileup,stderr=log)
	#success= 0
	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)
	open_mpileup.close()
	if not success:
		print "- Mpileup: %d min, %d sec" % elapsed_time
		return mpileup
	else:
		f.prRed('Error in Mpileup. Check log file.')
		exit(1)



#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     GATK     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def HaplotypeCaller(path,ram,bam,sample_name,reference,target,log,workdir):

	start_time = datetime.datetime.now()
	gvcf = workdir + '/'+ sample_name + '.g.vcf'
	args = ['java','-Xmx'+ram,'-jar',path,'-T','HaplotypeCaller','-R',reference,'-I', bam,'-o',gvcf,'-ERC','GVCF','--doNotRunPhysicalPhasing']

	if target != None:
		args += ['-L',target]

	success = subprocess.call(args,stdout=log,stderr=log)
	#success = 0
	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)
	
	if not success:
		print "- HaplotypeCaller "  +sample_name + ": %d min, %d sec" % elapsed_time
		return gvcf
	else:
		f.prRed('Error in gvcf generation: '+ sample_name+'. Check log file.')
		exit(1)

def GenotypeGVCFs(path,ram,name,gvcf_list,reference,target,log,workdir):

	start_time = datetime.datetime.now()
	#print "- GenotypeGVCFs"
	vcf = workdir + '/' + name + '.GATK.vcf'
	args = ['java','-Xmx'+ram,'-jar',path,'-T','GenotypeGVCFs','-R',reference,'-V:VCF', gvcf_list,'-o',vcf]

	if target != None:
		args += ['-L',target]

	success = subprocess.call(args,stdout=log,stderr=log)
	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)
	
	if not success:
		print "- GenotypeGVCFs: %d min, %d sec" % elapsed_time
		return vcf
	else:
		f.prRed('Error in GenotypeGVCFs. Check log file.')
		exit(1)



def Mutect2(path,ram,gbam,sbam,gsample_name,ssample_name,reference,log,workdir):

	start_time = datetime.datetime.now()
	vcf = workdir + '/' + name + '.Mutect.vcf'
	
	if gbam == None:
		args = ['java','-Xmx'+ram,'-jar',path,'-T','MuTect2','-R',reference,
			'-I:tumor',sbam,'-tumor',ssample_name,'-o',vcf]
	else:
		args = ['java','-Xmx'+ram,'-jar',path,'-T','MuTect2','-R',reference,
			'-I:tumor',sbam, '-I:normal',gbam,
			'-tumor',ssample_name,'-normal',nsample_name,'-o',vcf]

	if target != None:
		args += ['-L',target]

	success = subprocess.call(args,stdout=log,stderr=log)
	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)
	
	if not success:
		print "- Mutect2: %d min, %d sec" % elapsed_time
		return vcf
	else:
		f.prRed('Error in Mutect2. Check log file.')
		exit(1)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     FREEBAYES     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def FreeBayes(path,name,bam,bam_list,reference,target,log,workdir):

	start_time = datetime.datetime.now()
	#print "- FreeBayes"
	vcf = workdir + '/' + name + '.FreeBayes.vcf'
	if bam == None:
		args = [path,'-f',reference,'-L', bam_list,'-v',vcf,
			'--pooled-discrete','--pooled-continuous','--genotype-qualities','--report-genotype-likelihood-max','--allele-balance-priors-off']	

	else:
		args = [path,'-f',reference,'-b', bam,'-v',vcf,
			'--pooled-discrete','--pooled-continuous','--genotype-qualities','--report-genotype-likelihood-max','--allele-balance-priors-off']

	if target != None:
		args += ['-t',target]

	success = subprocess.call(args,stdout=log,stderr=log)
	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)
	
	if not success:
		print "- FreeBayes: %d min, %d sec" % elapsed_time
		return vcf
	else:
		f.prRed('Error in FreeBayes calling. Check log file.')
		exit(1)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     VARSCAN     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def VarScan_mpileup2snp(path,ram,mpileup,sample_list,reference,target,log,workdir):

	start_time = datetime.datetime.now()
	#print "- VarScan mpileup2snp"
	vcf = workdir + '/' + '.'.join(mpileup.split('/')[-1].split('.')[:-1] + ['VarScan.snp.vcf'])
	open_vcf = open(vcf,'w')
	args = ['java','-Xmx'+ram,'-jar',path,'mpileup2snp',mpileup,'--vcf-sample-list',sample_list,'--output-vcf','1','--strand-filter 0']

	if target != None:
		args += ['-L',target]

	#if filters == '0':
		#args += ['--min-coverage','1','--min-var-freq','0.01','--min-reads2','1','--min-avg-qual','0','--p-value','0.95']

	success = subprocess.call(args,stdout=open_vcf,stderr=log)
	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)
	open_vcf.close()
	
	if not success:
		print "- VarScan mpileup2snp: %d min, %d sec" % elapsed_time
		return vcf
	else:
		f.prRed('Error in VarScan mpileup2snp. Check log file.')
		exit(1)

def VarScan_mpileup2indel(path,ram,mpileup,sample_list,reference,target,log,workdir):

	start_time = datetime.datetime.now()
	#print "- VarScan mpileup2indel"
	vcf = workdir + '/' + '.'.join(mpileup.split('/')[-1].split('.')[:-1] + ['VarScan.indel.vcf'])
	open_vcf = open(vcf,'w')
	args = ['java','-Xmx'+ram,'-jar',path,'mpileup2indel',mpileup,'--vcf-sample-list',sample_list,'--output-vcf','1','--strand-filter 0']

	if target != None:
		args += ['-L',target]

	#if filters == '0':
		#args += ['--min-coverage','1','--min-var-freq','0.01','--min-reads2','1','--min-avg-qual','0','--p-value','0.95']

	success = subprocess.call(args,stdout=open_vcf,stderr=log)
	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)
	open_vcf.close()
	
	if not success:
		print "- VarScan mpileup2indel: %d min, %d sec" % elapsed_time
		return vcf
	else:
		f.prRed('Error in VarScan mpileup2indel. Check log file.')
		exit(1)

def VarScan_somatic(path,ram,mpileup,reference,target,log,workdir):

	start_time = datetime.datetime.now()
	#print "- VarScan mpileup2indel"
	vcf = workdir + '/' + '.'.join(mpileup.split('/')[-1].split('.')[:-1] + ['VarScan'])
	snp = vcf+'.snp.vcf'
	indel = svcf+'.indel.vcf'
	args = ['java','-Xmx'+ram,'-jar',path,'somatic',mpileup,'--output-vcf','1','--strand-filter 0','--mpileup 1']

	if target != None:
		args += ['-L',target]

	#if filters == '0':
		#args += ['--min-coverage','1','--min-var-freq','0.01','--min-reads2','1','--min-avg-qual','0','--p-value','0.95']

	success = subprocess.call(args,stdout=log,stderr=log)
	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)
	
	if not success:
		print "- VarScan mpileup2somatic: %d min, %d sec" % elapsed_time
		return snp,indel
	else:
		f.prRed('Error in VarScan mpileup2somatic. Check log file.')
		exit(1)

def Concat_VarScan_vcf(snp,indel,log):

	concat = '.'.join(snp.split('.')[:-2] + ['merge.vcf'])
	open_concat = open(concat,'w')
	sort = '.'.join(concat.split('.')[:-1] + ['sort.vcf'])
	open_sort = open(sort,'w')

	args = ['bgzip',snp]
	status = subprocess.call('bgzip -f '+snp, shell=True)
	status = subprocess.call('tabix -f '+snp+'.gz', shell=True)

	status = subprocess.call('bgzip -f '+indel, shell=True)
	status = subprocess.call('tabix -f '+indel+'.gz', shell=True)
	
	args = ['vcf-concat',snp+'.gz',indel+'.gz']
	success1 = subprocess.call(args,stdout=open_concat,stderr=log)

	args = ['vcf-sort','-c',concat]
	success2 = subprocess.call(args,stdout=open_sort,stderr=log)

	if success1 or success2:
		f.prRed('Error in VarScan concat. Check log file.')
		exit(1)
	else:
		return sort

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     VARDICT     ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


def Vardict(path,path_script,threads,gbam,sbam,gsample_name,ssample_name,design,reference,log,workdir):

	start_time = datetime.datetime.now()
	#print "- VarDict"
	vcf =  workdir + '/' + name + '.Vardict.vcf'
	open_vcf = open(vcf,'w')

	if sbam == None:
		args = [path,'-G',reference,'-f 0.05','-N',gsample_name,'-b',gbam]
		if threads != "":
			args += ['-th',threads]
		if target != None:
			args += ['-z 1','-F 0','-c 1','-S 2','-E 3','-g 4',target]
		args += ['|',path_script+'/teststrandbias.R','|',path_script+'/var2vcf_valid.pl','-f 0.05','-N',gsample_name]

	else:
		args = [path,'-G',reference,'-f 0.01','-N',ssample_name,'-b','"'+sbam+'|'+gbam+'"']
		if threads != "":
			args += ['-th',threads]
		if target != None:
			args += ['-z 1','-F 0','-c 1','-S 2','-E 3','-g 4',target]
		args += ['|',path_script+'/testsomatic.R','|',path_script+'/var2vcf_somatic.pl','-f 0.01','-N','"'+ssample_name+'|'+gsample_name+'"']

	success = subprocess.call(args,stdout=open_vcf,stderr=log)
	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)
	open_vcf.close()
	if not success:
		print "- VarDict: %d min, %d sec" % elapsed_time
		return vcf
	else:
		f.prRed('Error in VarDict. Check log file.')
		exit(1)

###------------------------------------------------------------------------------------------FEATURES EXTRACTION--------------------------------------------------------------------------------

def merge_vcf(path,name,gatk,freebayes,varscan,log,workdir):

	merge_vcf = workdir + '/' + name + '.merge.vcf'
	args = ['python',path,'-g',gatk,'-f',freebayes,'-v',varscan,'-o',merge_vcf]
	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		return merge_vcf
	else:
		f.prRed('Error in Vcf merging. Check log file.')
		exit(1)


def features_extractor(path,outpath,gatk,freebayes,varscan,merge,features_list,gvcf_path,design,log,workdir):
	start_time = datetime.datetime.now()
	args = ['python',path,'--listaFeatures',features_list,'--gvcf_path',gvcf_path,'-o',outpath]
	tsvfile = outpath+'/tsv.list'

	if merge != None:
		args += ['--merge',merge]
	else:
		args += ['-g',gatk,'-f',freebayes,'-v',varscan]

	if design == "Amplicon":
		args += ['-a']

	success = subprocess.call(args,stderr=log)
	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)

	if not success:
		#print "- Estraction: %d min, %d sec" % elapsed_time
		return merge,tsvfile
	else:
		f.prRed('Error in features extraction. Check log file.')
		exit(1)


def features_extractor_somatic(path,outpath,mutect,vardict,varscan,gname,sname,features_list,gvcf_path,design,log,workdir):
	start_time = datetime.datetime.now()
	args = ['python',path,'--listaFeatures',features_list,'--gvcf_path',gvcf_path,'-o',outpath]
	tsvfile = outpath+'/tsv.list'

	args += ['-g',mutect,'-f',vardict,'-v',varscan]

	args += ['-n',gname,'-t',sname]

	if design == "Amplicon":
		args += ['-a']

	success = subprocess.call(args,stderr=log)
	elapsed_time = divmod((datetime.datetime.now() - start_time).total_seconds(),60)

	if not success:
		#print "- Estraction: %d min, %d sec" % elapsed_time
		return merge,tsvfile
	else:
		f.prRed('Error in features extraction. Check log file.')
		exit(1)

def iEVA(path,vcf,reference,bam_list,log,workdir):

	ieva_vcf = workdir + '/' + '.'.join(mpileup.split('/')[-1].split('.')[:-1] + ['iEVA.vcf'])
	
	args = ['python',path,'--input',vcf,'--reference',reference,'--list',bam_list,'--outfile',ieva_vcf]

	args+=["--SimpleRepeat","--SimpleRepeatLength","--PseudoNucleotidesComposition","--RepeatMasker",
		"--gcContent","--VariantClass","--StrandBiasReads","--UnMappedReads","--MappingQualityZero",
		"--NotPrimaryAlignment","--SupplementaryAlignment","--NotPairedReads","--NotProperPairedReads",
		"--AlignmentScore","--NumberTotalDupReads","--NumberReadDupRef","--NumberReadDupAlt","--DuplicateReference",
		"--DuplicateAlternate","--DeltaDuplicate","--iEvaDepth","--iAlleleDepth","--ReadRef","--ReadAlt",
		"--MeanRefQscore","--MeanAltQscore","--TotalDPUnfilter","--NumberClippedReadsRef","--NumberClippedReadsAlt",
		"--ClippedReadsRef","--ClippedReadsAlt"]

	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		return ieva_vcf
	else:
		f.prRed('Error in iEVA annotation. Check log file.')
		exit(1)


###----------------------------------------------------------------------------------------------ANNOTATION-------------------------------------------------------------------------------------

def Vep(path,fork,name,reference,transcript_list,vep_af,vep_plugins,vcf,log,workdir):
	vep_vcf = workdir + '/' + name + '.VEP.vcf'
	
	args = [path,'-i',vcf,'-o',vep_vcf,'--species','homo_sapiens','--assembly','GRCh37','--force_overwrite','--no_stats',
		'--cache','--offline','--fasta',reference,'--merged',
		'--sift','b','--polyphen','b','--numbers','--vcf_info_field','ANN',
		'--hgvs','--transcript_version','--symbol','--canonical',
		'--check_existing','--vcf','--fork',fork,'--use_given_ref']

	if vep_af != "":
		for v in vep_af.split(','):
			args += [v]
	if vep_plugins['list'] != '':	
		args += add_vep_plugins(vep_plugins)	

	if transcript_list == None:
		args += ['--pick']
	if transcript_list == "most_severe":
		args += ['--most_severe']
	
	#print args
	success = subprocess.call(args,stdout=log,stderr=log)

	if not success:
		print "- VEP"
		return vep_vcf
	else:
		f.prRed('Error in VEP annotation. Check log file.')
		exit(1)

def add_vep_plugins(vep_plugins):
	args = []
	plugins = vep_plugins['list'].split(',')
	for plugin in plugins:
		pl_args = plugin
		if vep_plugins[plugin]['path'] != "":
			pl_args += ',' + vep_plugins[plugin]['path']
		if vep_plugins[plugin]['fields'] != "":	
			pl_args += ',' + vep_plugins[plugin]['fields']

		args += ['--plugin',pl_args]
	return args


def add_Annotation(path,name,vcf,tsv,tag_list,transcript_list,log,workdir):

	annotated_tsv = workdir + '/' + name + '.tsv'
	
	args = ['python',path,'--vcf',vcf,'--tsv',tsv,'-o',annotated_tsv,'--tag_list',tag_list,'--trs_list',transcript_list]
	success = subprocess.call(args,stdout=log,stderr=log)

	if success:
		f.prRed('Error in tsv annotation. Check log file.')
		exit(1)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::     TOOLS     :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

def header_fix(path,vcf,variantcaller,log):
	fix_vcf = '.'.join(vcf.split('.')[:-1] + ['hfix.vcf'])
	open_fix_vcf = open(fix_vcf,'w')
	args = ['python',path,'-v',variantcaller,'-f',vcf]
	success = subprocess.call(args,stdout=open_fix_vcf,stderr=log)
	open_fix_vcf.close()
	if not success:
		return fix_vcf
	else:
		f.prRed('Error in Vcf fix. Check log file.')
		exit(1)

def vcf_norm(path,vcf,reference,log):
	norm_vcf = '.'.join(vcf.split('.')[:-1] + ['norm.vcf'])
	args = [path,'norm','-m','-both','-d','both','-f',reference,vcf,'-o',norm_vcf]
	success = subprocess.call(args,stdout=log,stderr=log)
	
	if not success:
		return norm_vcf
	else:
		f.prRed('Error in Vcf normalization. Check log file.')
		exit(1)


