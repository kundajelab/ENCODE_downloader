#!/usr/bin/env python3
import sys, os, time
import json
import requests
import subprocess
import collections
import re
import argparse

parser = argparse.ArgumentParser(prog='ENCODE fastq downloader', \
                                    description='Download fastqs from the ENCODE portal and \
                                    generate Kundaje lab BDS pipeline shell script for all samples. \
                                    If authenticaion information (--encode-access-key-id and --encode-secret-key) is given, \
                                    unpublished fastqs only visible to submitters with valid authentication \
                                    can be downloaded.')
parser.add_argument('dir_download', metavar='dir-download', type=str, \
                        help='Root directory to save downloaded fastqs, directory structure: \
                            dir_download/award_rfa/assay_title/assay_category/each_sample_accession_id. \
                            ')
parser.add_argument('award_rfa', metavar='award-rfa', type=str, \
                        help='Award RFA (e.g. ENCODE3)')
parser.add_argument('assay_category', metavar='assay-category', type=str, \
                        help='Assay category (e.g. DNA+accessibility)')
parser.add_argument('assay_title', metavar='assay-title', type=str, \
                        help='Assay title (e.g. ATAC-seq)')
parser.add_argument('scientific_name', metavar='scientific-name', type=str, \
                        help='Scientific name for genome (e.g. Mus+musculus, Homo+sapiens)')
group_accession_ids = parser.add_mutually_exclusive_group()
group_accession_ids.add_argument('--ignored-accession-ids-file', type=str, \
                        help='Accession IDs in this text file will be ignored. (1 acc. ID per line)')
group_accession_ids.add_argument('--accession-ids-file', type=str, \
                        help='Only accession IDs in this text file will be downloaded. (1 acc. ID per line). Others will be ignored.')
parser.add_argument('--max-download', type=int, default=8, \
                        help='Maximum number of fastqs for concurrent downloading. Forced to 1 if used with --encode-access-key-id')
parser.add_argument('--encode-access-key-id', type=str, \
                        help='To see files only visible to submitters vith a valid authentication. \
                        Get your access key ID and secret key from the portal homepage menu YourID->Profile->Add Access key')
parser.add_argument('--encode-secret-key', type=str, \
                        help='ENCODE secret key (--encode-access-key-id must be specified).' )
parser.add_argument('--dry-run', action="store_true",\
                        help='The downloader shows a list of fastqs to be downloaded but does not download them.')
parser.add_argument('--encode-url-base', type=str, default='https://www.encodeproject.org', \
                        help='URL base for the ENCODE portal (e.g. https://www.encodeproject.org')
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
# parser.add_argument('--pipeline-extra-parameters', type=str, default='', \
#                         help='Extra parameters will be appended to the BDS pipeline command line. Use \\- for -')

args = parser.parse_args()

# commonly used string
mid_underscore = args.award_rfa+'_'+args.assay_category+'_'+args.assay_title+'_'+args.scientific_name
mid_slash = args.award_rfa+'/'+args.assay_category+'/'+args.assay_title+'/'+args.scientific_name

# loaded ignored accession list

ignored_accession_ids = []
if args.ignored_accession_ids_file and os.path.isfile(args.ignored_accession_ids_file):
    with open(args.ignored_accession_ids_file,'r') as f:
        ignored_accession_ids = f.read().splitlines()
    ignored_accession_ids = \
        [accession_id for accession_id in ignored_accession_ids if accession_id and not accession_id.startswith("#") ]
    print '* ignored_accession_ids:\n', ignored_accession_ids

accession_ids = []
if args.accession_ids_file and os.path.isfile(args.accession_ids_file):
    with open(args.accession_ids_file,'r') as f:
        accession_ids = f.read().splitlines()
    accession_ids = \
        [accession_id for accession_id in accession_ids if accession_id and not accession_id.startswith("#") ]
    print '* accession_ids:\n', accession_ids

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
HEADERS = {'accept': 'application/json'}
encode_search_url = args.encode_url_base+'/search/?format=json' \
                    +'&type=Experiment' \
                    +'&limit=all' \
                    +'&assay_slims='+args.assay_category \
                    +'&assay_title='+args.assay_title \
                    +'&award.rfa='+args.award_rfa \
                    +'&replicates.library.biosample.donor.organism.scientific_name='+args.scientific_name
print('* query:', encode_search_url)

if args.encode_access_key_id and not args.encode_secret_key or \
    not args.encode_access_key_id and args.encode_secret_key:
    print("both --encode-access-key-id and --encode-secret-key must be specified.")    
    raise ValueError

if args.encode_access_key_id: # if ENCODE key is given
    encode_auth = (args.encode_access_key_id, args.encode_secret_key)
    search_data = requests.get(encode_search_url, headers=HEADERS, auth=encode_auth)
else:
    search_data = requests.get(encode_search_url, headers=HEADERS)    

# print search_data
json_data_search = search_data.json() #json.loads(search_data)     
cnt_accession = 0
for item in json_data_search['@graph']:
    accession_id = item['accession']
    url_suffix = item['@id']
    # get json from ENCODE portal for accession id
    if args.encode_access_key_id: # if ENCODE key is given
        search_data = requests.get(args.encode_url_base+url_suffix+'?format=json',headers=HEADERS, auth=encode_auth)
    else:
        search_data = requests.get(args.encode_url_base+url_suffix+'?format=json',headers=HEADERS)
    json_data_exp = search_data.json()
    # if accession_id == "ENCSR229QKB" : print(json_data_exp)
    print('* %s' % (accession_id,))
    if ignored_accession_ids and accession_id in ignored_accession_ids: # ignore if in the black list
        print('ignored')
        continue
    elif accession_ids and not accession_id in accession_ids:
        print('not in the list')
        continue
    cnt_accession += 1
    # for pipeline script
    cmd_pipeline = '# %d\nTITLE=%s; WORKDIR=%s; mkdir -p $WORKDIR; cd $WORKDIR\n' \
                        % (cnt_accession,accession_id,args.pipeline_run_dir+'/'+mid_slash+'/'+accession_id)
    param_pipeline = ''
    
    fastqs = dict()
    for org_f in json_data_exp['original_files']:
        if args.encode_access_key_id: # if ENCODE key is given
            search_fastq = requests.get(args.encode_url_base+org_f+'?format=json',headers=HEADERS,auth=encode_auth)
        else:
            search_fastq = requests.get(args.encode_url_base+org_f+'?format=json',headers=HEADERS)
        f = search_fastq.json()
        if f['status']=='error': continue
        pair = int(f['paired_end']) if 'paired_end' in f else -1
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
        if not args.dry_run:
            os.system(cmd_mkdir)

        # check number of downloading fastqs
        while int(subprocess.check_output('ps aux | grep wget | wc -l', shell=True).strip('\n')) \
                    > args.max_download-1:
            print '# of downloads exceeded the limit (%d), retrying in 20 seconds...' % (args.max_download,)
            time.sleep(20)

        # download fastq
        if args.encode_access_key_id:
            basename = url_fastq.split("/")[-1]
            filename = '%s/%s' % (dir,basename)
            if os.path.exists(filename):
                print('File exists: %s, rep:%d, pair:%d' % (url_fastq, bio_rep_id, pair))
            else:
                cmd_curl = 'curl -RL -u %s:%s %s -o %s' % (args.encode_access_key_id, \
                        args.encode_secret_key, url_fastq, filename)
                if args.dry_run:
                    print('Dry-run: %s, rep:%d, pair:%d' % (url_fastq, bio_rep_id, pair))
                else:
                    print('Downloading: %s, rep:%d, pair:%d' % (url_fastq, bio_rep_id, pair))
                    os.system(cmd_curl)
                # print(cmd_curl)
            # else:
                # print("already exists")
        else:
            cmd_wget = 'wget -bqcN -P %s %s' % (dir,url_fastq)
            if args.dry_run:
                print('Dry-run: %s, rep:%d, pair:%d' % (url_fastq, bio_rep_id, pair))
            else:
                print('Downloading: %s, rep:%d, pair:%d' % (url_fastq, bio_rep_id, pair))
                os.system(cmd_wget)
                time.sleep(0.25) # wait for 0.25 second per fastq

        # check if paired with other fastq                
        paired_with = None
        if 'paired_with' in f:
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
    # if args.pipeline_extra_parameters: param_ENCODE_meta += args.pipeline_extra_parameters

    param_award_rfa = '-' + args.award_rfa + ' ' # -ENCODE3

    cmd_pipeline += 'bds_scr $TITLE %s %s %s %s %s \n' \
                    % (args.pipeline_script, param_basic, param_pipeline, param_ENCODE_meta, param_award_rfa )                    
    cmd_pipeline += 'sleep 0.5\n\n'
    #print cmd_pipeline
    file_pipeline.write(cmd_pipeline)
