#!/usr/bin/env python2
import sys, os, time
import json
import urllib2
import subprocess
import collections
import re
import argparse

parser = argparse.ArgumentParser(prog='ENCODE fastq downloader', \
                                    description='Download fastqs from the ENCODE portal and generate Kundaje lab BDS pipeline shell script for all samples.')
parser.add_argument('dir_download', metavar='dir-download', type=str, \
                        help='Root directory to save downloaded fastqs')
parser.add_argument('award_rfa', metavar='award-rfa', type=str, \
                        help='Award RFA (e.g. ENCODE3)')
parser.add_argument('assay_category', metavar='assay-category', type=str, \
                        help='Assay category (e.g. DNA+accessibility)')
parser.add_argument('assay_title', metavar='assay-title', type=str, \
                        help='Assay title (e.g. ATAC-seq)')
parser.add_argument('scientific_name', metavar='scientific-name', type=str, \
                        help='Scientific name for genome (e.g. Mus+musculus, Homo+sapiens)')
parser.add_argument('--encode-url-base', type=str, default='https://www.encodeproject.org', \
                        help='URL base for the ENCODE portal (e.g. https://www.encodeproject.org')
parser.add_argument('--ignored-accession-ids-file', type=str, \
                        help='Accession IDs in this text file will be ignored. (1 acc. ID per line)')
parser.add_argument('--max-download', type=int, default=8, \
                        help='Maximum number of fastqs for concurrent downloading')
parser.add_argument('--pipeline-num-thread', type=str, default='${PIPELINE_NTH}', \
                        help='Number of threads for each pipeline')
parser.add_argument('--pipeline-run-dir', type=str, default='${PIPELINE_RUN_DIR}', \
                        help='Root directory of outputs for pipeline')
parser.add_argument('--pipeline-data-dir', type=str, default='${PIPELINE_DATA_DIR}', \
                        help='Root directory of data for pipeline (recommended to match with the first parameter [DIR_DOWNLOAD])')
parser.add_argument('--pipeline-script', type=str, default='${PIPELINE_SCRIPT}', \
                        help='BDS pipeline script (e.g. /path/to/atac.bds)')
parser.add_argument('--pipeline-web-url-base', type=str, default='${PIPELINE_WEB_URL_BASE}', \
                        help='URL base for browser tracks (e.g. http://mitra.stanford.edu/kundaje')
parser.add_argument('--pipeline-ref-genome', type=str, default='${REF_GENOME}', \
                        help='Reference genome name for pipeline (e.g. hg38_ENCODE3, mm10_ENCODE3, hg19, mm9)')
parser.add_argument('--pipeline-encode-lab', type=str, default='', \
                        help='ENCODE lab for pipeline (e.g. /labs/anshul-kundaje/)')
parser.add_argument('--pipeline-encode-award', type=str, default='', \
                        help='ENCODE award for pipeline (e.g. /awards/U41HG007000/)')
parser.add_argument('--pipeline-encode-assembly', type=str, default='', \
                        help='ENCODE assembly (ref. genome name in ENCODE database) for pipeline (e.g. hg38 (x), GRCh38 (o)')
parser.add_argument('--pipeline-encode-alias-prefix', type=str, default='', \
                        help='ENCODE alias prefix for pipeline (pipeline output files will have aliases of [prefix].[filename], lab name is recommended, e.g. anshul-kundaje)')

args = parser.parse_args()

# commonly used string
mid_underscore = args.award_rfa+'_'+args.assay_category+'_'+args.assay_title+'_'+args.scientific_name
mid_slash = args.award_rfa+'/'+args.assay_category+'/'+args.assay_title+'/'+args.scientific_name

# loaded ignored accession list
ignored_accession_ids = []
if args.ignored_accession_ids_file and os.path.isfile(args.ignored_accession_ids_file):
    ignored_accession_ids = open(args.ignored_accession_ids_file,'r').read().splitlines()
print '* ignored_accession_ids:'
print [accession_id for accession_id in ignored_accession_ids if accession_id and not accession_id.startswith("#") ]

# init shell script for pipeline
# if os.path.exists(args.pipeline_run_dir):
#     file_pipeline = open(args.pipeline_run_dir+'/run_pipelines_'+mid_underscore+'.bash','w')
# elif os.path.exists(args.dir_download):
#     file_pipeline = open(args.dir_download+'/run_pipelines_'+mid_underscore+'.sh','w')
# else:
#     file_pipeline = open('./run_pipelines_'+mid_underscore+'.sh','w')
file_pipeline = open('./run_pipelines_'+mid_underscore+'.sh','w')

file_pipeline.write('#!/bin/bash\n')
file_pipeline.write('# specify your own pipeline run/data directories\n')
file_pipeline.write('pipeline_run_dir='+args.pipeline_run_dir+'\n')
file_pipeline.write('DIR_PIPELINE_DATA='+args.pipeline_data_dir+'\n')
file_pipeline.write('# path for pipeline script (atac.bds, chipseq.bds ...)\n')
file_pipeline.write('BDS_PIPELINE_SCRIPT='+args.pipeline_script+'\n')
file_pipeline.write('# URL base for genome browser tracks\n')
file_pipeline.write('WEB_URL_BASE='+args.pipeline_web_url_base+'\n\n')
    
# send query to ENCODE portal and parse
encode_search_url = args.encode_url_base+'/search/?format=json' \
                    +'&type=Experiment' \
                    +'&limit=all' \
                    +'&assay_slims='+args.assay_category \
                    +'&assay_title='+args.assay_title \
                    +'&award.rfa='+args.award_rfa \
                    +'&replicates.library.biosample.donor.organism.scientific_name='+args.scientific_name
print '* query:', encode_search_url

try:
    # request = urllib2.Request(encode_search_url)
    # request.add_header('User-Agent','Mozilla/5.0')
    # search_data = urllib2.build_opener().open(request).read()
    search_data = urllib2.urlopen(encode_search_url).read()
except urllib2.HTTPError, e:
    print e.code
    print e.msg
    exit()

# search_data = subprocess.check_output('curl -H "Accept: application/json" "%s"' \
#                         % (encode_search_url,), shell=True ).strip('\n')
# print search_data
json_data_search = json.loads(search_data)     
# print json_data_search
cnt_accession = 0
for item in json_data_search['@graph']:
    accession_id = item['accession']
    url_suffix = item['@id']
    # get json from ENCODE portal for accession id
    json_data_exp=json.loads(urllib2.urlopen(args.encode_url_base+url_suffix+'?format=json').read())
    print '* %s' % (accession_id,)
    if accession_id in ignored_accession_ids: # ignore if in the black list
        print 'ignored'
        continue
    cnt_accession += 1
    # for pipeline script
    cmd_pipeline = '# %d\nTITLE=%s; WORKDIR=%s; mkdir -p $WORKDIR; cd $WORKDIR\n' \
                        % (cnt_accession,accession_id,args.pipeline_run_dir+'/'+mid_slash+'/'+accession_id)
    param_pipeline = ''
    
    fastqs = dict()
    for f in json_data_exp['files']:
        pair = int(f['paired_end']) if f.has_key('paired_end') else -1
        if 'fastq' != f['file_type']: # ignore if not fastq
            continue
        url_fastq = args.encode_url_base+f['href']
        bio_rep_id = int(f['replicate']['biological_replicate_number'])
        tech_rep_id = int(f['replicate']['technical_replicate_number'])

        if tech_rep_id > 1:
            break;
        # create directory for downloading
        dir_suffix = '/'+mid_slash+'/'+accession_id+'/rep'+str(bio_rep_id)
        if pair > 0:
            dir_suffix += '/pair'+str(pair)
        dir = args.dir_download + dir_suffix
        cmd_mkdir = 'mkdir -p %s' % (dir,)
        os.system(cmd_mkdir)

        # check number of downloading fastqs
        while int(subprocess.check_output('ps aux | grep wget | wc -l', shell=True).strip('\n')) \
                    > args.max_download-1:
            print '# of downloads exceeded the limit (%d), retrying in 20 seconds...' % (args.max_download,)
            time.sleep(20)

        # download fastq
        print 'Downloading: %s, rep:%d, pair:%d' % (url_fastq, bio_rep_id, pair)
        cmd_wget = 'wget -bqcN -P %s %s' % (dir,url_fastq)
        os.system(cmd_wget)
        # wait for 0.25 second per fastq
        time.sleep(0.25)

        # check if paired with other fastq                
        paired_with = None
        if f.has_key('paired_with'):
            paired_with = f['paired_with'].split('/')[2]

        # relative path for fastq (for pipeline)
        rel_fastq = args.pipeline_data_dir + dir_suffix +'/'+os.path.basename(url_fastq)
        accession_id_fastq = os.path.basename(url_fastq).replace('.fastq.gz','')                

        # store fastqs with the same bio_rep_id and pair: these fastqs will be pooled later in a pipeline                                
        fastqs[accession_id_fastq] = (bio_rep_id, pair, paired_with, rel_fastq)

    cnt_fastq_to_be_pooled = collections.defaultdict(int)
    already_done = []
    for accession_id_fastq in fastqs:
        bio_rep_id, pair, paired_with, rel_fastq = fastqs[accession_id_fastq]
        if rel_fastq in already_done:
            continue                
        cnt_fastq_to_be_pooled[bio_rep_id] += 1
        # suffix for fastqs to be pooled
        param_endedness = ('-pe' if paired_with else '-se' ) + str(bio_rep_id)
        if cnt_fastq_to_be_pooled[bio_rep_id] > 1:
            suffix = '_'+str(cnt_fastq_to_be_pooled[bio_rep_id])
            suffix2 = ':'+str(cnt_fastq_to_be_pooled[bio_rep_id])
            param_endedness = ''
        else:
            suffix = ''
            suffix2 = ''
        if paired_with:
            _, pair2, _, rel_fastq2 = fastqs[paired_with]
            cmd_pipeline += 'FASTQ%d_%d%s=%s\n' % (bio_rep_id, pair, suffix, rel_fastq)
            cmd_pipeline += 'FASTQ%d_%d%s=%s\n' % (bio_rep_id, pair2, suffix, rel_fastq2)
            param_pipeline += ' %s -fastq%d_%d%s $FASTQ%d_%d%s -fastq%d_%d%s $FASTQ%d_%d%s' \
                                % ( param_endedness, bio_rep_id, pair, suffix2, bio_rep_id, pair, suffix, \
                                               bio_rep_id, pair2, suffix2, bio_rep_id, pair2, suffix )
            already_done.append(rel_fastq2)
        else:
            cmd_pipeline += 'FASTQ%d%s=%s\n' % (bio_rep_id, suffix, rel_fastq)
            param_pipeline += ' %s -fastq%d%s $FASTQ%d%s' \
                                % (param_endedness, bio_rep_id, suffix2, bio_rep_id, suffix)
        already_done.append(rel_fastq)

    param_basic = '-nth %s -title $TITLE -species %s -url_base %s ' \
                    % (args.pipeline_num_thread, args.pipeline_ref_genome, \
                        args.pipeline_web_url_base+'/'+mid_slash+'/'+accession_id+'/out' )
    param_ENCODE_meta = '-ENCODE_accession %s -ENCODE_assay_category %s -ENCODE_assay_title %s -ENCODE_award_rfa %s ' \
                    % (accession_id, args.assay_category, args.assay_title, args.award_rfa )
    if args.pipeline_encode_lab: param_ENCODE_meta += '-ENCODE_lab '+args.pipeline_encode_lab+' '
    if args.pipeline_encode_award: param_ENCODE_meta += '-ENCODE_award '+args.pipeline_encode_award+' '
    if args.pipeline_encode_assembly: param_ENCODE_meta += '-ENCODE_assembly '+args.pipeline_encode_assembly+' '
    if args.pipeline_encode_alias_prefix: param_ENCODE_meta += '-ENCODE_alias_prefix '+args.pipeline_encode_alias_prefix+' '

    param_award_rfa = '-' + args.award_rfa + ' ' # -ENCODE3

    cmd_pipeline += 'bds_scr $TITLE %s %s %s %s %s \n' \
                    % (args.pipeline_script, param_basic, param_pipeline, param_ENCODE_meta, param_award_rfa )                    
    cmd_pipeline += 'sleep 0.5\n\n'
    #print cmd_pipeline
    file_pipeline.write(cmd_pipeline)
