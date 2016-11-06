#!/usr/bin/env python2
import sys, os, time
import json
import urllib2
import subprocess
import collections
import re
import argparse

parser = argparse.ArgumentParser(prog='ENCODE fastq downloader', \
                                    description='Download fastqs from the ENCODE portal.', \
                                    epilog='python2 ENCODE_downloader.py '+ \
                                            '/path/download ENCODE3 DNA+accessibility '+ \
                                            'ATAC-seq Homo-sapiens --ref-genome hg38_ENCODE3')
parser.add_argument('dir_download', metavar='dir-download', type=str, \
                        help='Root directory to save downloaded fastqs')
parser.add_argument('award_rfa', metavar='award-rfa', type=str, \
                        help='Award RFA (e.g. ENCODE3)')
parser.add_argument('assay_category', metavar='assay-category', type=str, \
                        help='Assay category (e.g. DNA+accessibility)')
parser.add_argument('assay_title', metavar='assay-title', type=str, \
                        help='Assay title (e.g. ATAC-seq)')
parser.add_argument('scientific_name', metavar='scientific-name', type=str, \
                        help='Scientific name for genome (e.g. Mus+musculus, Homo-sapiens)')
parser.add_argument('--ignored-accession-ids-file', type=str, \
                        help='Accession IDs in this text file will be ignored. (1 acc. ID per line)')
parser.add_argument('--max-download', type=int, default=8, \
                        help='Maximum number of fastqs for concurrent downloading')
parser.add_argument('--ref-genome', type=str, default='${REF_GENOME}', \
                        help='Reference genome name for pipeline (e.g. hg38_ENCODE3, mm10_ENCODE3, hg19, mm9)')
parser.add_argument('--num-thread-pipeline', type=int, default='${NTH_PIPELINE}', \
                        help='Number of thread for each pipeline')
parser.add_argument('--dir-pipeline-run', type=str, default='${DIR_PIPELINE_RUN}', \
                        help='Root directory of outputs for pipeline')
parser.add_argument('--dir-pipeline-data', type=str, default='${DIR_PIPELINE_DATA}', \
                        help='Root directory of data for pipeline')
parser.add_argument('--bds-pipeline-script', type=str, default='${BDS_PIPELINE_SCRIPT}', \
                        help='BDS pipeline script (e.g. /path/to/atac.bds)')
parser.add_argument('--web-url-base', type=str, default='${WEB_URL_BASE}', \
                        help='URL base for browser tracks (e.g. http://mitra.stanford.edu/kundaje')
parser.add_argument('--encode-url-base', type=str, default='https://www.encodeproject.org', \
                        help='URL base for the ENCODE portal (e.g. https://www.encodeproject.org')
args = parser.parse_args()

# commonly used string
mid_underscore = args.award_rfa+'_'+args.assay_category+'_'+args.assay_title+'_'+args.scientific_name
mid_slash = args.award_rfa+'/'+args.assay_category+'/'+args.assay_title+'/'+args.scientific_name

# init shell script for pipeline
if os.path.exists(args.dir_pipeline_run):
    file_pipeline = open(args.dir_pipeline_run+'/run_pipelines_'+mid_underscore+'.sh','w')
elif os.path.exists(args.dir_download):
    file_pipeline = open(args.dir_download+'/run_pipelines_'+mid_underscore+'.sh','w')
else:
    file_pipeline = open('./run_pipelines_'+mid_underscore+'.sh','w')
file_pipeline.write('# specify your own pipeline run/data directories\n')
file_pipeline.write('DIR_PIPELINE_RUN='+args.dir_pipeline_run+'\n')
file_pipeline.write('DIR_PIPELINE_DATA='+args.dir_pipeline_data+'\n')
file_pipeline.write('# path for pipeline script (atac.bds, chipseq.bds ...)\n')
file_pipeline.write('BDS_PIPELINE_SCRIPT='+args.bds_pipeline_script+'\n')
file_pipeline.write('# URL base for genome browser tracks\n')
file_pipeline.write('WEB_URL_BASE='+args.web_url_base+'\n\n')

# loaded ignored accession list
ignored_accession_ids = []
if os.path.isfile(args.ignored_accession_ids_file):
    ignored_accession_ids = open(args.ignored_accession_ids_file,'r').read().splitlines()
print 'ignored_accession_ids:\n', ignored_accession_ids

# send query to ENCODE portal and parse
encode_search_url = args.encode_url_base+'/search/?type=Experiment' \
                    +'&limit=all' \
                    +'&assay_slims='+args.assay_category \
                    +'&assay_title='+args.assay_title \
                    +'&award.rfa='+args.award_rfa \
                    +'&replicates.library.biosample.donor.organism.scientific_name='+args.scientific_name
print(encode_search_url)
search_data = urllib2.urlopen(encode_search_url).read()
for line in search_data.split('\n'):
    if line.startswith('{'):
        json_data_search = json.loads(line)     
        cnt_accession = 0
        for item in json_data_search['@graph']:
            accession_id = item['accession']
            url_suffix = item['@id']
            # get json from ENCODE portal for accession id
            json_data_exp=json.loads(urllib2.urlopen(args.encode_url_base+url_suffix+'?format=json').read())
            print '===== %s =====' % (accession_id,)
            if accession_id in ignored_accession_ids: # ignore if in the black list
                print 'ignored'
                continue
            cnt_accession += 1
            # for pipeline script
            cmd_pipeline = '# %d\nTITLE=%s; WORKDIR=%s; mkdir -p $WORKDIR; cd $WORKDIR\n' \
                                % (cnt_accession,accession_id,args.dir_pipeline_run+'/'+mid_slash+'/'+accession_id)
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
                #os.system(cmd_wget)
                # wait for 0.25 second per fastq
                #time.sleep(0.25)

                # check if paired with other fastq                
                paired_with = None
                if f.has_key('paired_with'):
                    paired_with = f['paired_with'].split('/')[2]

                # relative path for fastq (for pipeline)
                rel_fastq = args.dir_pipeline_data + dir_suffix +'/'+os.path.basename(url_fastq)
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
            
            param_award = '-' + args.award_rfa # -ENCODE3
            cmd_pipeline += 'bds_scr $TITLE %s -nth %d %s -url_base %s -title $TITLE -species %s %s\n' \
                            % (args.bds_pipeline_script, args.num_thread_pipeline, param_award, args.web_url_base, \
                                args.ref_genome, param_pipeline)
            cmd_pipeline += 'sleep 1\n\n'
            #print cmd_pipeline
            file_pipeline.write(cmd_pipeline)

"""
python2 ENCODE_downloader.py /srv/scratch/shared/surya/leepc12/data ENCODE3 DNA+accessibility ATAC-seq Homo+sapiens \
--ref-genome hg38_ENCODE3 --bds-pipeline-script '${PIPELINE_CODE}/atac_dnase_pipelines/atac.bds' --num-thread-pipeline 3 \
--web-url-base http://mitra.stanford.edu/kundaje --ignored-accession-ids-file ignored_accession_ids.txt
"""
