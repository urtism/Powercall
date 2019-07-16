import argparse
import re

def estrai_annotazione():
    vcf = open(opts.vcf,'r')
    out = open(opts.out,'w')
    tsv = open(opts.tsv,'r')
    tags = open(opts.tag_list,'r')
    var_list = dict()
    tag_list=[]
    
    for var in tsv:
        var = var.rstrip()

        if var.startswith('CHROM'):
            header=var

        else:
            chrom = var.split('\t')[0]
            pos = var.split('\t')[1]
            id = var.split('\t')[2]
            ref = var.split('\t')[3]
            alt =  var.split('\t')[4]
            
            var_id = '\t'.join([chrom,pos,ref,alt])
            var_list[var_id] = var.split('\t')

    for tag in tags:
        tag=tag.rstrip()
        if tag.startswith('#'):
            continue
        else:
            tag_list+=[tag]

    tags.close()
    header_tot = header.split('\t') + tag_list
    out.write('\t'.join(header_tot)+'\n')

    for line in vcf:
        line = line.rstrip()

        if line.startswith('##INFO=<ID=ANN') or line.startswith('##INFO=<ID=CSQ'):
            start = line.find('Allele')
            end = line.find('">')
            header_ann = (line[start:end]).split('|')
            continue
        elif line.startswith('#'):
            continue
            
        else:
            chrom = line.split('\t')[0]
            pos = line.split('\t')[1]
            id = line.split('\t')[2]
            ref = line.split('\t')[3]
            alt =  line.split('\t')[4]
            var_id = '\t'.join([chrom,pos,ref,alt])
            
            if var_id in var_list.keys():
                info = (line.split('\t')[7]).split('=')[1]
                info_split = info.split('|')

                tags = var_list.get(var_id)
                
                for tag in tag_list:              
                    new_tag = info_split[header_ann.index(tag)]               
                    if new_tag == '':
                        new_tag = '.'
                    tags = tags + [new_tag]
                
                var_list[var_id] = tags
                out.write('\t'.join(var_list.get(var_id))+'\n')


def tags_extractor(tag_list):
    tag_l = open(tag_list,'r')
    tags = []
    for tag in tag_l:
        if not tag.startswith('#'):
            tags += [tag.rstrip()]
    tag_l.close()
    return tags

def transcr_extractor(transcr_list):
    transcr_l = open(transcr_list,'r')
    transcrs = []
    for transcr in transcr_l:
        
        if not transcr.startswith('#'):
            transcrs += [transcr.rstrip()]
    transcr_l.close()
    return transcrs

def split_annotation(anninfo,header,transcrs):
    ann_array = [[],[],[]]
    for ann in anninfo.split(','):
        ann_split = ann.split('|')
        tr = ann_split[header.index('Feature')].split('.')[0]
        if tr in transcrs:
            #print tr
            ann_array[0] = ann_split

        if ann_split[header.index('CANONICAL')]:
            ann_array[1] = ann_split
        
        ann_array[2] += [ann_split]

    if ann_array[0] == []:
        ann_array[0] = ann_array[2][0]

    return ann_array
    


def main():
    parser = argparse.ArgumentParser('aggiunge le annotazioni fornite nel file -f (una per riga) al file di input -i in un formato tab delimited in output -o')
    
    parser.add_argument('-i','--vcf',help="file delle varianti annotate in formato vcf")
    parser.add_argument('-l','--tag_list',help="lista di annotation features da aggiungere: una per riga",default=None)
    parser.add_argument('-f','--tsv',help="file tab delimited a cui aggiungere le annotation features",default=None)
    parser.add_argument('-t','--trs_list',help="lista di trascritti su cui splittare le annotazioni",default=None)
    parser.add_argument('-o','--out',help="file di output")
    parser.add_argument('-O','--other_transcripts',help="file di output containing variants in other transcripts",default=None)
    parser.add_argument('-S', '--split', help="abilita lo split delle varianti per i trascritti contenuti in list", action='store_true')
    
    global opts
    
    opts = parser.parse_args()
    
    
    vcf = open(opts.vcf,'r')
    out = open(opts.out,'w')
    tsv = open(opts.tsv,'r')

    if opts.other_transcripts:
        other = open(opts.other_transcripts,'w')

    tag_list = tags_extractor(opts.tag_list)
    print 'trs list ' +opts.trs_list 
    if opts.trs_list != None and opts.trs_list != '':
        transcrs = transcr_extractor(opts.trs_list)
    else:
        transcrs = []

    var_list = dict()
    
    for var in tsv:
        var = var.rstrip()

        if var.startswith('CHROM'):
            header=var

        else:
            chrom = var.split('\t')[0]
            pos = var.split('\t')[1]
            id = var.split('\t')[2]
            ref = var.split('\t')[3]
            alt =  var.split('\t')[4]
            
            var_id = '\t'.join([chrom,pos,ref,alt])
            var_list[var_id] = var.split('\t')

    for line in vcf:
        line = line.rstrip()

        if line.startswith('##INFO=<ID=ANN') or line.startswith('##INFO=<ID=CSQ'):
            start = line.find('Allele')
            end = line.find('">')
            header_ann = (line[start:end]).split('|')
            header_tot = header.split('\t') + tag_list
            out.write('\t'.join(header_tot)+'\n')
            try:
                other.write('\t'.join(header_tot)+'\n')
            except:
                continue
            continue
        elif line.startswith('#'):
            continue
        else:
            tags = []
            chrom = line.split('\t')[0]
            pos = line.split('\t')[1]
            id = line.split('\t')[2]
            ref = line.split('\t')[3]
            alt =  line.split('\t')[4]
            var_id = '\t'.join([chrom,pos,ref,alt])
            info = line.split('\t')[7]
            for i in info.split(';'):
                if i.startswith('ANN=') or i.startswith('CSQ='):
                    anninfo = re.sub('ANN=','',i)

            if var_id in var_list.keys():

                annotations = split_annotation(anninfo,header_ann,transcrs)
                
                if opts.trs_list != None:
                    main_trs = annotations[0]
                    other_trs = annotations[1] + annotations[2]
                else:
                    main_trs = annotations[1]
                
                for tag in tag_list:
                    try:           
                        new_tag = main_trs[header_ann.index(tag)]
                    except:
                        print var_id,annotations
                    if new_tag == '':
                        new_tag = '.'
                    tags = tags + [new_tag]
                
                toprint = var_list[var_id] + tags
                out.write('\t'.join(toprint)+'\n')

                if opts.other_transcripts != None:
                    for trs in other_trs:
                        tags = []
                        for tag in tag_list:
                            new_tag = main_trs[header_ann.index(tag)]

                            if new_tag == '':
                                new_tag = '.'
                            tags = tags + [new_tag]
                        
                        toprint = var_list[var_id] + tags
                        other.write('\t'.join(toprint)+'\n')
    vcf.close()     
    out.close()
    tsv.close()
main()