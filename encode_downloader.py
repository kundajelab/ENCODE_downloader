#!/usr/bin/env python
'''
Written by Jin Lee (leepc12@gmail.com) 
Anshul Kundaje's lab at Stanford University.
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
                            help='List of file types to be downloaded: fastq (default). \
                            You can define file_type, file_type:file_format:output_type or \
                            file_type:output_type. \
                            Supported file_type: fastq, bam, bed, bigWig, bigBed. \
                            Supported file_format: fastq, bam, bed, bigWig, bigBed. \
                            Output file_format: alignments, unfiltered alignments, ... \
                            For details, checkout ENCODE portal\'s schema. \
                            e.g. --file-types fastq "bam:unfiltered alignments" "bam:alignments".')
    parser.add_argument('--assemblies', nargs='+', default=['all'], type=str,
                            help='Assemblies (reference used for mapping) allowed. e.g. --assemblies GRCh38 hg19.')
    parser.add_argument('--encode-access-key-id', type=str,
                            help='To download unpublished files visible to submitters vith a valid authentication. \
                            Get your access key ID and secret key from the portal homepage menu YourID->Profile->Add Access key')
    parser.add_argument('--encode-secret-key', type=str,
                            help='ENCODE secret key (--encode-access-key-id must be specified).' )
    parser.add_argument('--ignored-accession-ids-file', type=str,
                            help='Text file with ignored accession IDs.')    
    parser.add_argument('--pooled-rep-only', action="store_true",\
                            help='Download genome data from pooled replicates only.')
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
    # convert hg38->GRCh38
    for i, assembly in enumerate(args.assemblies):
        if assembly == 'hg38':
            args.assemblies[i] = 'GRCh38'
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
                if debug: print('{}: {}'.format(help_str, str(json)))
                ret |= True
        except UnicodeEncodeError:
            pass
    return ret

def get_depth_one( json_obj ):
    result = {}
       # add info to metadata json
    for key in json_obj:
        val = json_obj[key]
        if not type(val)==dict and not type(val)==list:
            result[key] = val
    return result

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
    print(accession_ids)

    os.system('mkdir -p {}'.format(args.dir))
    # ordered dict to write metadata table (including all accessions)
    all_file_metadata = collections.OrderedDict()
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
                    # print('Warning: inferred assembly ({}) from keys and values in json (e.g. organism name)...'.format(assembly))
                    break
        if 'assay_category' in json_data_exp:        
            assay_category = json_data_exp['assay_category']
        else:
            assay_category = None
        # read files in accession
        downloaded_this_exp_accession = False
        
        # init metadata object
        metadata = get_depth_one(json_data_exp)
        metadata['files'] = {} # file info
        for org_f in json_data_exp['original_files']:
            # if args.encode_access_key_id: # if ENCODE key is given
            #     search_file = requests.get(ENCODE_BASE_URL+org_f+'?format=json',headers=HEADERS,auth=encode_auth)
            # else:
            #     search_file = requests.get(ENCODE_BASE_URL+org_f+'?format=json',headers=HEADERS)
            retry_cnt = 0
            while True:
                try:
                    if args.encode_access_key_id: # if ENCODE key is given
                        search_file = requests.get(ENCODE_BASE_URL+org_f+'?format=json',headers=HEADERS,auth=encode_auth)
                    else:
                        search_file = requests.get(ENCODE_BASE_URL+org_f+'?format=json',headers=HEADERS)
                except:
                    print('Exception caught, retrying in 120 seconds...')
                else:
                    break
                retry_cnt += 1
                if retry_cnt>100:
                    raise Exception('Exceeded maximum number of retries {}. Aborting...'.format(retry_cnt-1))
                print('Retrial: {}'.format(retry_cnt))
                time.sleep(120)

            f = search_file.json()
            status = f['status'].lower().replace(' ','_')
            if status=='error': continue
            file_assembly = f['assembly'] if 'assembly' in f else ''
            arr = f['file_type'].lower().split(' ')
            if len(arr)>1:
                file_type = arr[0]
                file_format = arr[1]
            else:
                file_type = arr[0]
                file_format = file_type                
            output_type = f['output_type'].lower() #.replace(' ','_')
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
            if type(tech_rep_id)==list and len(tech_rep_id)>0:
                tech_rep_id=tech_rep_id[0]
            else:
                tech_rep_id=tech_rep_id            
            if args.pooled_rep_only and type(bio_rep_id)==list and len(bio_rep_id)<2:
                continue
            # print(file_accession_id, file_assembly, file_type, file_format, output_type, bio_rep_id, tech_rep_id, pair)
            if file_type == 'fastq':
                # if tech_rep_id != '1' and tech_rep_id != 1: break;
                pass
            else:
                if not 'all' in args.assemblies and not file_assembly in args.assemblies: continue
            # create directory for downloading
            dir_suffix = accession_id+'/'+status+'/'+file_assembly+'/'+output_type.replace(' ', '_')+'/'+file_type.replace(' ', '_')
            if file_type!=file_format: dir_suffix += '/'+file_format
            # if bio_rep_id>0: dir_suffix += '/rep'+str(bio_rep_id)
            if bio_rep_id:
                dir_suffix += '/rep'+'_rep'.join([str(i) for i in bio_rep_id])
            if pair>0: dir_suffix += '/pair'+str(pair)
            dir = args.dir+'/' + dir_suffix

            # print file_accession_id, file_assembly, file_type, file_format, output_type, bio_rep_id, tech_rep_id, pair, dir
            # continue

            if not args.dry_run:
                os.system('mkdir -p {}'.format(dir))
            # continue
            # check number of downloading files
            # while int(subprocess.check_output('ps aux | grep wget | wc -l', shell=True).strip('\n')) \
            while int(subprocess.check_output('ps aux | grep wget | wc -l', shell=True)) \
                        > args.max_download-1:
                print('# of downloads exceeded the limit ({}), retrying in 20 seconds...'.format(args.max_download))
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
            rel_file = rel_file.replace('//','/')

            # check if paired with other fastq                
            paired_with = None
            if file_type == 'fastq' and 'paired_with' in f:
                paired_with = f['paired_with'].split('/')[2]

            # for fastq, store files with the same bio_rep_id and pair: these files will be pooled later in a pipeline
            if bio_rep_id:                
                metadata['files'][file_accession_id] = dict(
                    file_type=file_type,
                    file_format=file_format,
                    output_type=output_type,
                    status=status,
                    bio_rep_id=bio_rep_id,
                    pair=pair,
                    paired_with=paired_with,
                    rel_file=rel_file)
            downloaded_this_exp_accession = True

        if not args.dry_run and downloaded_this_exp_accession:
            os.system('mkdir -p {}'.format(args.dir+'/'+accession_id))
            with open(args.dir+'/'+accession_id+'/metadata.org.json',mode='w') as fp:
                fp.write(json.dumps(json_data_exp, indent=4))
            with open(args.dir+'/'+accession_id+'/metadata.json',mode='w') as fp:
                fp.write(json.dumps(metadata, indent=4))
            all_file_metadata[accession_id] = metadata['files']

    # make TSV for all downloaded files
    if not args.dry_run and all_file_metadata:
        # count max. number of files per exp. accession
        max_num_files = max( [len(all_file_metadata[acc_id]) for acc_id in all_file_metadata] )
        # table header
        header = 'accession\t'+\
            'description(comma-delimited; file_acc_id:status:file_type:file_format:output_type:bio_rep_id:pair,...)\tfile'+ \
            '\tfile'.join([str(i+1) for i in range(max_num_files)]) + '\n'
        contents = ''
        # table contents        
        for accession_id in all_file_metadata:
            file_metadata = all_file_metadata[accession_id]
            desc = ''
            tmp_cnt = 0
            for file_acc_id in file_metadata:
                tmp_cnt += 1
                metadata = file_metadata[file_acc_id]
                desc += ':'.join(   [file_acc_id,
                                    metadata['status'],
                                    metadata['file_type'],
                                    metadata['file_format'],
                                    metadata['output_type'],
                                    metadata['file_type'],
                                    '_'.join(str(x) for x in metadata['bio_rep_id']),
                                    str(metadata['pair']) ])
                if tmp_cnt<len(file_metadata):
                    desc += ','
            files = '\t'.join( [file_metadata[a]['rel_file'] for a in file_metadata] )
            contents += '\t'.join([accession_id,desc,files]) + '\n'
        with open(args.dir+'/all_files.tsv',mode='w') as fp:
            fp.write(header)
            fp.write(contents)

if __name__=='__main__':
    main()
