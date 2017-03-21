#!/usr/bin/env python
'''
    code written by Jin Lee (leepc12@gmail.com) in Anshul Kundaje's lab at Stanford University.
'''

import sys
import os
import time
import json
import requests
import subprocess
import collections
import re
import argparse
import ENCODE_kundaje_pipeline

ENCODE_BASE_URL = 'https://www.encodeproject.org'

def parse_arguments():
    parser = argparse.ArgumentParser(prog='ENCODE downloader',
                                        description='Downloads genome data files from the ENCODE portal \
                                        and save them under [WORK_DIR]/[ACCESSION_ID]/. \
                                        Also generates Kundaje lab BDS pipeline shell script for all samples (for fastq and bam only). \
                                        If authentication information (--encode-access-key-id and --encode-secret-key) is given, \
                                        unpublished files only visible to submitters with valid authentication \
                                        can be downloaded.')
    parser.add_argument('url_or_file', metavar='url-or-file', nargs='+', type=str,
                            help='List of URLs/files/accession_ids \
                                (ENCODE search/experiment URL, exp. accesion ids text file or exp. accession id). \
                                Make sure that URL is quotted. \
                                e.g. ENCSR000ELE exp_acc_ids.txt "https://www.encodeproject.org/search/?\
                                type=Experiment&assay_term_name=DNase-seq\
                                &replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens\
                                &biosample_term_name=CD14-positive+monocyte\
                                &month_released=October%%2C+2011".\
                            ')
    parser.add_argument('--dir', default='.', type=str,
                            help='[WORK_DIR] : Root directory for all downloaded genome data.')
    parser.add_argument('--file-types', nargs='+', default=['fastq'], type=str,
                            help='List of file types to be downloaded: fastq (default), bam, bed, bigWig, bigBed, ... \
                            e.g. --file-types fastq bam.')
    parser.add_argument('--assemblies', nargs='+', default=['all'], type=str,
                            help='Assemblies (reference used for mapping) allowed. e.g. --assemblies GRCh38 hg19.')
    parser.add_argument('--encode-access-key-id', type=str,
                            help='To download unpublished files visible to submitters vith a valid authentication. \
                            Get your access key ID and secret key from the portal homepage menu YourID->Profile->Add Access key')
    parser.add_argument('--encode-secret-key', type=str,
                            help='ENCODE secret key (--encode-access-key-id must be specified).' )
    parser.add_argument('--ignored-accession-ids-file', type=str,
                            help='Text file with ignored accession IDs.')    
    parser.add_argument('--dry-run', action="store_true",\
                            help='Dry-run: downloads nothing, but generates pipeline shell script.')
    parser.add_argument('--dry-run-list-accession-ids', action="store_true",
                            help='Dry-run: downloads nothing, but show a list of accession IDs matching URL.')
    parser.add_argument('--max-download', type=int, default=8,
                            help='Maximum number of files for concurrent downloading. \
                            Parallel downloading will be disabled when used with --encode-access-key-id')
    parser.add_argument('--assembly-map', nargs='+', default=['Mus+musculus:mm10','Homo+sapiens:GRCh38'], type=str,
                            help='List of strings to infer ENCODE assembly from species name; [SPECIES_NAME]:[ASSEMBLY]. \
                            e.g. --assembly-map Mus+musculus:mm10 Homo+sapiens:GRCh38')
    parser.add_argument('--pipeline-atac-bds-path', type=str, default='${ATAC_BDS_DIR}/atac.bds',
                            help='Path for atac.bds.')
    parser.add_argument('--pipeline-chipseq-bds-path', type=str, default='${CHIPSEQ_BDS_DIR}/chipseq.bds',
                            help='Path for chipseq.bds.')
    parser.add_argument('--pipeline-web-url-base', type=str,
                            help='URL base for browser tracks. e.g. http://mitra.stanford.edu/kundaje')
    parser.add_argument('--pipeline-encode-lab', type=str, default='',
                            help='ENCODE lab for pipeline. e.g. /labs/anshul-kundaje/')
    parser.add_argument('--pipeline-encode-award', type=str, default='',
                            help='ENCODE award for pipeline. e.g. /awards/U41HG007000/')
    parser.add_argument('--pipeline-encode-award-rfa', type=str, default='',
                            help='ENCODE award RFA for pipeline. e.g. ENCODE3')
    parser.add_argument('--pipeline-encode-alias-prefix', type=str, default='',
                            help='ENCODE alias prefix for pipeline (pipeline output files will have aliases of [prefix]:[filename], \
                            lab name is recommended). e.g. anshul-kundaje)')
    group_ignore_status = parser.add_mutually_exclusive_group()
    group_ignore_status.add_argument('--ignore-released', action='store_true', \
                            help='Ignore released data (except fastqs).')
    group_ignore_status.add_argument('--ignore-unpublished', action='store_true', \
                            help='Ignore unpublished data.')
    args = parser.parse_args()

    if args.encode_access_key_id and not args.encode_secret_key or \
        not args.encode_access_key_id and args.encode_secret_key:
        print("Both parameters --encode-access-key-id and --encode-secret-key must be specified together.")
        raise ValueError
    args.dir = os.path.abspath(args.dir)
    # make file_types lowercase
    for i, file_type in enumerate(args.file_types):
        args.file_types[i] = file_type.lower()
    return args

def is_encode_url( url ):
    return url.startswith(ENCODE_BASE_URL)

def is_encode_search_query_url( url ):
    return url.startswith(ENCODE_BASE_URL+'/search/?')

def is_encode_exp_url( url ):
    return url.startswith(ENCODE_BASE_URL+'/experiments/ENCSR')

def get_accession_id_from_encode_exp_url( url ):    
    for s in url.split('/')[-2:]:
        if s.startswith('ENC'): return s
    return None

def is_file_type_fastq( file_type ):
    return file_type.lower() in ['fastq']

def is_file_type_bam( file_type ):
    return file_type.lower() in ['bam']

def get_accession_ids( accession_ids_file ):
    accession_ids = []
    if accession_ids_file and os.path.isfile(accession_ids_file):
        with open(accession_ids_file,'r') as f:
            accession_ids = f.read().splitlines()
        accession_ids = [accession_id for accession_id in accession_ids \
                if accession_id and accession_id.strip() \
                    and not accession_id.strip().startswith("#") \
                    and not accession_id.strip().startswith("/") ]
    return accession_ids

def deep_search(json, s, help_str='root',debug=False):
    ret = False
    if type(json)==list:
        for i, e in enumerate(json):
            ret |= deep_search(e, s, help_str+'[{}]'.format(i))
    elif type(json)==dict:
        for k in json:        
            ret |= deep_search(json[k], s, help_str+'.{}'.format(k))
    else:
        try:
            if s in str(json):
                if debug: print '{}: {}'.format(help_str, str(json))
                ret |= True
        except UnicodeEncodeError:
            pass
    return ret

def main():
    args = parse_arguments()

    # read ignored accession ids
    ignored_accession_ids = get_accession_ids( args.ignored_accession_ids_file )

    HEADERS = {'accept': 'application/json'}
    if args.encode_access_key_id: # if ENCODE key is given
        encode_auth = (args.encode_access_key_id, args.encode_secret_key)
    accession_ids = []
    # process multiple inputs
    for url_or_file in args.url_or_file:
        if is_encode_search_query_url(url_or_file):
            if not 'limit=all' in url_or_file:
                url_or_file += '&limit=all'
            if not 'format=json' in url_or_file:
                url_or_file += '&format=json'
            # send query to ENCODE portal and parse
            if args.encode_access_key_id: # if ENCODE key is given
                search_data = requests.get(url_or_file, headers=HEADERS, auth=encode_auth)
            else:
                search_data = requests.get(url_or_file, headers=HEADERS)    
            json_data_search = search_data.json() #json.loads(search_data)
            for item in json_data_search['@graph']:
                accession_id = item['accession']
                accession_ids.append(accession_id)
        elif is_encode_exp_url(url_or_file):
            accession_id = get_accession_id_from_encode_exp_url(url_or_file)
            accession_ids.append( accession_id )
        elif os.path.exists(url_or_file) and os.path.isfile(url_or_file):
            accession_ids += get_accession_ids( url_or_file )
        elif url_or_file.startswith('ENCSR'):
            accession_ids.append(url_or_file)
        else:
            print("Only URL, accession_ids_file or accession_id is allowed for input ({}).".format(url_or_file))
            raise ValueError
    print accession_ids

    pipeline_sh = ENCODE_kundaje_pipeline.PipelineShellScript( args.dir, args.pipeline_encode_lab, 
        args.pipeline_encode_alias_prefix, args.pipeline_encode_award, os.path.abspath(args.dir), 
        args.pipeline_web_url_base, args.pipeline_chipseq_bds_path, args.pipeline_atac_bds_path)

    os.system('mkdir -p {}'.format(args.dir))

    # download files for each accession id
    for accession_id in accession_ids:
        # get accession info
        print("="*10+" "+accession_id+" "+"="*10)
        if args.dry_run_list_accession_ids: continue
        if ignored_accession_ids and accession_id in ignored_accession_ids: # ignore if in the black list
            print('\tignored (--ignored-accession-ids-file)')
            continue
        if accession_ids and not accession_id in accession_ids:
            print('\tnot in the list (--accession-ids-file)')
            continue
        # get json from ENCODE portal for accession id
        if args.encode_access_key_id: # if ENCODE key is given
            search_data = requests.get(ENCODE_BASE_URL+'/experiments/'+accession_id+'?format=json',
                headers=HEADERS, auth=encode_auth)
        else:
            search_data = requests.get(ENCODE_BASE_URL+'/experiments/'+accession_id+'?format=json',
                headers=HEADERS)
        json_data_exp = search_data.json()  

        if json_data_exp['status']=='error':
            print("Error: cannot access to accession {}".format(accession_id))
            print(json_data_exp)
            continue
        assembly = None
        if not assembly: # infer assembly from organism name...
            for s in args.assembly_map:
                arr = s.replace('+',' ').split(':')
                if deep_search(json_data_exp, arr[0]):
                    assembly = arr[1]
                    print('Warning: inferred assembly ({}) from keys and values in json (e.g. organism name)...'.format(assembly))
                    break

        assay_title = json_data_exp['assay_title']
        if 'assay_category' in json_data_exp:        
            assay_category = json_data_exp['assay_category']
        else:
            assay_category = None
        award_rfa = args.pipeline_encode_award_rfa
        # write metadata file on each accession dir.
        # with open(args.dir+'/'+accession_id+'/metadata.txt',mode='w') as f: 
        #     f.write('assay_category={}\n'.format(assay_category))
        #     f.write('assay_title={}\n'.format(assay_title))
        #     f.write('assembly={}\n'.format(assembly))
        dir = args.dir+'/'+accession_id
        cmd_mkdir = 'mkdir -p {}'.format(dir)
        if not args.dry_run:
            os.system(cmd_mkdir)
        with open(dir+'/metadata.json',mode='w') as f:
            f.write(json.dumps(json_data_exp, indent=4))
        # read files in accession
        files = {}
        for org_f in json_data_exp['original_files']:
            if args.encode_access_key_id: # if ENCODE key is given
                search_file = requests.get(ENCODE_BASE_URL+org_f+'?format=json',headers=HEADERS,auth=encode_auth)
            else:
                search_file = requests.get(ENCODE_BASE_URL+org_f+'?format=json',headers=HEADERS)
            f = search_file.json()
            status = f['status'].lower().replace(' ','_')
            if status=='error': continue
            arr = f['file_type'].lower().split(' ')
            if len(arr)>1:
                file_type = arr[0]
                file_format = arr[1]
            else:
                file_type = arr[0]
                file_format = file_type                
            output_type = f['output_type'].lower().replace(' ','_')
            valid = ('all' in args.file_types)
            for ft in args.file_types:
                arr = ft.split(':')
                if len(arr)>2:
                    if file_type==arr[0] and file_format==arr[1] and output_type==arr[2]:
                        valid = True
                        break                    
                elif len(arr)>1:
                    if (file_type==arr[0] or file_format==arr[0]) and \
                        (file_format==arr[1] or output_type==arr[1]):
                        valid = True
                        break
                else:                    
                    if file_type==arr[0]:
                        valid = True
                        break
            if not valid: continue            
            url_file = ENCODE_BASE_URL+f['href']
            file_accession_id = f['accession']

            if args.ignore_released and status=='released': continue
            if args.ignore_unpublished and status!='released': continue
            if 'paired_end' in f:
                pair = int(f['paired_end']) 
            else:
                pair = -1
            if 'replicate' in f and 'biological_replicates_number'in f['replicate']:
                bio_rep_id = f['replicate']['biological_replicate_number']
            else: # 'biological_replicates'in f:
                bio_rep_id = f['biological_replicates']
            if 'replicate' in f and 'technical_replicate_number'in f['replicate']:
                tech_rep_id = f['replicate']['technical_replicate_number']
            else: #if 'technical_replicates' in f:
                tech_rep_id = f['technical_replicates']
            if len(bio_rep_id)>1:
                bio_rep_id=0
            else:
                bio_rep_id=int(bio_rep_id[0])
            if type(tech_rep_id)==list:
                tech_rep_id=tech_rep_id[0]
            else:
                tech_rep_id=tech_rep_id
            file_assembly = f['assembly'] if 'assembly'in f else ''
            # print file_accession_id, file_assembly, file_type, file_format, output_type, bio_rep_id, tech_rep_id, pair, dir

            if file_type == 'fastq':
                # if tech_rep_id != '1' and tech_rep_id != 1: break;
                pass
            else:
                if not 'all' in args.assemblies and not file_assembly in args.assemblies: break;
            # create directory for downloading
            dir_suffix = accession_id+'/'+status+'/'+file_assembly+'/'+output_type+'/'+file_type
            if file_type!=file_format: dir_suffix += '/'+file_format
            if bio_rep_id>0: dir_suffix += '/rep'+str(bio_rep_id)            
            if pair>0: dir_suffix += '/pair'+str(pair)
            dir = args.dir+'/' + dir_suffix

            # print file_accession_id, file_assembly, file_type, file_format, output_type, bio_rep_id, tech_rep_id, pair, dir
            # continue

            cmd_mkdir = 'mkdir -p {}'.format(dir)
            if not args.dry_run:
                os.system(cmd_mkdir)
            # continue
            # check number of downloading files
            while int(subprocess.check_output('ps aux | grep wget | wc -l', shell=True).strip('\n')) \
                        > args.max_download-1:
                print '# of downloads exceeded the limit ({}), retrying in 20 seconds...'.format(args.max_download)
                time.sleep(20)

            # download file
            basename = url_file.split("/")[-1]
            filename = '{}/{}'.format(dir,basename)
            if os.path.exists(filename):
                print('File exists ({}): {}, rep:{}, pair:{}'.format(file_type, url_file, bio_rep_id, pair))
            elif args.dry_run:
                print('Dry-run ({}): {}, rep:{}, pair:{}'.format(file_type, url_file, bio_rep_id, pair))
            else:
                print('Downloading ({}): {}, rep:{}, pair:{}'.format(file_type, url_file, bio_rep_id, pair))
                if args.encode_access_key_id:
                    cmd_curl = 'curl -RL -u {}:{} {} -o {}'.format(args.encode_access_key_id,
                            args.encode_secret_key, url_file, filename)
                    os.system(cmd_curl)
                else:
                    cmd_wget = 'wget -bqcN -P {} {}'.format(dir,url_file)
                    os.system(cmd_wget)
                    time.sleep(0.25) # wait for 0.25 second per fastq

            # relative path for file (for pipeline)            
            rel_file = args.dir + '/' + dir_suffix + '/' + os.path.basename(url_file)

            # check if paired with other fastq                
            paired_with = None
            if file_type == 'fastq' and 'paired_with' in f:
                paired_with = f['paired_with'].split('/')[2]

            # for fastq, store files with the same bio_rep_id and pair: these files will be pooled later in a pipeline
            files[file_accession_id] = (file_type, status, bio_rep_id, pair, paired_with, rel_file)

        # print accession_id, assembly, assay_title, award_rfa, files
        pipeline_sh.add_sample(accession_id, assembly, assay_title, assay_category, award_rfa, files)

    pipeline_sh.write_to_file()

if __name__=='__main__':
    main()