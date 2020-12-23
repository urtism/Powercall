

import re
import string
import argparse
import sys
import statistics
import scipy.stats as stats

def calcola_gt(GT):

	if '0/1' in GT or '1/0' in GT or '0|1' in GT or '1|0' in GT:
		return '0/1'

	elif '1/1' in GT:
		return '1/1'

	else:
		return '0/0'


def calcola_ad(AD):
	alt_c = []
	ref_c = []
	for ad in AD:
		if ad != '.':
			ref_c += [int(ad.split(',')[0])]
			alt_c += [int(ad.split(',')[1])]

	try:
		ref = stats.median(ref_c)
		alt = stats.median(alt_c)
	except:
		ref = '.'
		alt = '.'

	return ','.join([str(ref),str(alt)])


def calcola_dp(DP):

	try:
		DP.remove('.')
	except:
		print(DP)
		return '.'
	try:	
		return stats.median(DP)
	except:
		print(DP)
		return '.'


def calcola_af(AF):

	try:
		AF.remove('.')
	except:
		print(AF)
		return '.'
	try:	
		return stats.mean(AF)
	except:
		print(AF)
		return '.'


if __name__ == '__main__':
	
	parser = argparse.ArgumentParser('Adapt VCF to eVAI, extracting information about input sample')
	parser.add_argument('-v', '--vcf', help="input vcf")
	parser.add_argument('-o', '--outdir', help="output ")
	parser.add_argument('-s', '--sample', help="sample to be extract, output vcf will be named as outdir/sample.vcf")

	
	global opts 
	opts = parser.parse_args()

	vcf = open(opts.vcf,'r').readlines()
	out_vcf = open(opts.outdir + '/' opts.sample +'.vcf'+,'w')
	
	new_header = ['##fileformat=VCFv4.2',
		'##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
		'##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fractions of alternate alleles in the tumor">',
		'##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
		'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">']

	sample_name = opts.sample

	new_header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+ sample_name
	out_vcf.write('\n'.join(new_header)+'\n')
	
	for line in vcf:
		if line.startswith('##'):
			continue
		elif line.startswith('#CHROM'):
			header = line.rstrip().split('\t')
		else:
			var = line.rstrip().split('\t')
			chrom, pos, id, ref, alt, qual, filter, info, format = line.rstrip().split('\t')[:8]

			sample_format = var[header.index(sample_name)].split(':')

			GT = []
			AD = []
			DP = []
			AF = []

			for field in format.split(':'):
				if field.startswith('GT_'):
					GT += [sample_format[format.split(':').index(field)]]
					gt = calcola_gt(GT)

				elif field.startswith('AD_'):

					AD += [sample_format[format.split(':').index(field)]]

				elif field.startswith('DP_'):

					DP += [sample_format[format.split(':').index(field)]]

				elif field.startswith('AF_'):

					AF += [sample_format[format.split(':').index(field)]]


			gt = calcola_gt(GT)
			ad = calcola_ad(AD)
			dp = calcola_dp(DP)
			af = calcola_af(AF)

			new_format_sample = ':'.join([gt, ad, dp, af])

			qual = '.'
			filter = '.'
			info = '.'
			format = 'GT:AD:AF:DP'

		out_vcf.write('\t'.join([chrom, pos, id, ref, alt, qual, filter, info, format, new_format_sample])+'\n')
