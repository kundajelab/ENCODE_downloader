#!/usr/bin/env python

import json
import os
import sys
import argparse
import collections

PIPELINE_SH_ITEM_TEMPLATE = '''
# SN={sn}
TITLE={title}
WORKDIR={pipeline_out_root_dir}/{title}; mkdir -p $WORKDIR; cd $WORKDIR
bds_scr {title} {pipeline_bds_script} -title {title} -species {species} \\
{input_file_param}
{input_end_param}
'''

def parse_arguments():
    parser = argparse.ArgumentParser(prog='Kundaje lab pipeline BDS shell script generator',
                        description='Generates run_pipelines.sh')
    parser.add_argument('--exp-acc-ids-file', type=str,
                            help='File with experiment accession id in each line.')
    parser.add_argument('--exp-data-root-dir', type=str,
                            help='Root directory where you downloaded experiment data.')
    parser.add_argument('--ctl-data-root-dir', type=str,
                            help='Root directory where you downloaded control data.')
    parser.add_argument('--exp-id-to-ctl-id-file', type=str,
                            help='File with exp_id[TAB]ctl_id in each line. Only for samples with controls.')
    parser.add_argument('--pipeline-sh-filename', type=str, default="run_pipelines.sh",
                            help='Path for atac.bds or chipseq.bds.')
    parser.add_argument('--pipeline-out-root-dir', type=str, default='', 
                            help='Pipeline output root directory.')
    parser.add_argument('--pipeline-bds-script', type=str, required=True,
                            help='Path for atac.bds or chipseq.bds.')
    parser.add_argument('--exp-file-type', type=str, required=True,
                            choices=['fastq','bam','filt_bam'],
                            help='Choose file_type of experiment replicates to run pipelines. \
                                Files with other types will be ignored exp. replicates.')
    parser.add_argument('--ctl-file-type', type=str,
                            choices=['fastq','bam','filt_bam'],
                            help='Optinoal. If not defined then --exp-file-type will be useds for controls. \
                                Choose file_type of control replicates to run pipelines. \
                                Files with other types will be ignored for controls.')
    args = parser.parse_args()

    if args.ctl_data_root_dir and not args.exp_id_to_ctl_id_file or \
        not args.ctl_data_root_dir and args.exp_id_to_ctl_id_file:
        raise Exception('--ctl-data-root-dir and --exp-id-to-ctl-id-file must be defined together.')
    args.pipeline_out_root_dir = os.path.abspath(args.pipeline_out_root_dir)
    ctl_exists = args.ctl_data_root_dir!=None
    return args, ctl_exists

def mkdir_p(path):
    if os.path.exists(path): return
    os.makedirs(path)

def read_acc_ids_file(f):
    acc_ids=[]
    with open(f,'r') as fp:
        lines = fp.readlines()
        for line in lines:
            acc_ids.append(line.strip())
    return acc_ids

def read_exp_to_ctl_file(f):
    map_exp_to_ctl={}
    with open(f,'r') as fp:
        lines = fp.readlines()
        for line in lines:
            arr = line.strip().split('\t')
            exp_id = arr[0]
            ctl_id = arr[1]
            map_exp_to_ctl[exp_id] = ctl_id
    return map_exp_to_ctl

def match_file_type(file_type, output_type, file_type_to_run_pipeline):
    if file_type_to_run_pipeline=='fastq' and file_type!='fastq':
        return False
    if file_type_to_run_pipeline=='bam' and \
        (file_type!='bam' or output_type!='unfiltered alignments') :
        return False
    if file_type_to_run_pipeline=='filt_bam' and \
        (file_type!='bam' or output_type!='alignments') :
        return False
    return True

def parse_file_acc_json(obj):
    file_type = obj['file_type']
    output_type = obj['output_type']
    bio_rep_id = obj['bio_rep_id'][0]
    paired_with = obj['paired_with']
    pair = obj['pair']
    rel_file = obj['rel_file']
    return file_type, output_type, bio_rep_id, pair, paired_with, rel_file

def is_paired_end(metadata_org_json_file):
    with open(metadata_org_json_file,'r') as fp:
        json_obj = json.load(fp)
    for f_obj in json_obj['files']:
        if 'run_type' in f_obj:
            return f_obj['run_type']=='paired-ended'
    raise Exception('could not find endedness information from {}'.format(
        metadata_org_json_file))

def infer_species(metadata_org_json_file):
    with open(metadata_org_json_file,'r') as fp:
        json_obj = json.load(fp)
    for assembly in json_obj['assembly']:
        if assembly=='GRCh38' or assembly=='hg38': return 'hg38'
        if assembly=='hg19': return 'hg19'
        if assembly=='GRCm38' or assembly=='mm10': return 'mm10'
        if assembly=='mm9': return 'mm9'
    if deep_search(json_obj, 'Homo sapiens'):
        return 'hg38'
    if deep_search(json_obj, 'Mus Musculus'):
        return 'mm10'
    raise Exception('could not find/infer species from {}'.format(
        metadata_org_json_file))

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

def parse_metadata_json_file(json_file, file_type_to_run_pipeline):
    with open(json_file,'r') as fp:
        json_obj = json.load(fp)
    result = []
    files = json_obj['files']

    # bio_rep_id does not always start from rep1 sometimes it's like [rep3, rep5]
    # so make it start from rep1 and increment [rep3, rep5] -> [rep1, rep2]
    bio_rep_ids = collections.defaultdict(int)
    map_bio_rep_id_to_serial_bio_rep_id = {}
    
    bio_rep_id_pairs = collections.defaultdict(int)
    fastq_merge_id = {}
    for file_acc_id in files:
        file_type, output_type, bio_rep_id, pair, paired_with, rel_file \
            = parse_file_acc_json(files[file_acc_id])
        if not match_file_type(file_type, output_type, file_type_to_run_pipeline):
            continue
        if file_type=='fastq' and pair==2:
            continue
        bio_rep_ids[bio_rep_id] += 1

        bio_rep_id_pairs[(bio_rep_id,pair)] += 1
        fastq_merge_id[file_acc_id] = bio_rep_id_pairs[(bio_rep_id,pair)]

        print(file_type, output_type, bio_rep_id, pair, paired_with, rel_file)

    # sort bio_rep_ids
    cnt = 0
    for bio_rep_id in sorted(bio_rep_ids):
        cnt += 1
        map_bio_rep_id_to_serial_bio_rep_id[bio_rep_id] = cnt

    # iterate over all file accession ids
    # map new serial bio_rep_id    
    for file_acc_id in files:
        file_type, output_type, bio_rep_id, pair, paired_with, rel_file \
            = parse_file_acc_json(files[file_acc_id])
        if not match_file_type(file_type, output_type, file_type_to_run_pipeline):
            continue
        if file_type=='fastq' and pair==2:
            continue

        if paired_with: # fastq pe
            merge_id = fastq_merge_id[file_acc_id]

            file_type2, output_type2, bio_rep_id2, pair2, paired_with2, rel_file2 \
                = parse_file_acc_json(files[paired_with])
            (file_acc_id,paired_with)
            new_bio_rep_id2 = map_bio_rep_id_to_serial_bio_rep_id[bio_rep_id2]
            result.append((file_type2, output_type2,
                new_bio_rep_id2, pair2, merge_id, paired_with2, rel_file2))
        else: # fastq se or other file_type
            merge_id = 1

        new_bio_rep_id = map_bio_rep_id_to_serial_bio_rep_id[bio_rep_id]
        result.append((file_type, output_type, 
            new_bio_rep_id, pair, merge_id, paired_with, rel_file))
    return result

def parse_exp_metadata_json(exp, ctl):
    input_file_param = ''

    for (file_type, output_type, bio_rep_id, pair, merge_id,
        paired_with, rel_file) in exp:

        if file_type=='fastq':
            input_file_param += '-fastq{}_{}{} {} \\\n'.format(
                bio_rep_id,
                pair if paired_with else 1,
                ':{}'.format(merge_id) if merge_id>1 else '',
                rel_file)
        elif file_type=='bam' \
            and output_type=='alignments':
            input_file_param += '-filt_bam{} {} \\\n'.format(
                bio_rep_id,                
                rel_file)

        elif file_type=='bam' \
            and output_type=='unfiltered alignments':
            input_file_param += '-bam{} {} \\\n'.format(
                bio_rep_id,                
                rel_file)
        else:
            Exception('fastq and bam input only!')

    if ctl:
        for (file_type, output_type, bio_rep_id, merge_id,
            pair, paired_with, rel_file) in ctl:

            if file_type=='fastq':
                input_file_param += '-ctl_fastq{}_{}{} {} \\\n'.format(
                    bio_rep_id,
                    pair if paired_with else 1,
                    ':{}'.format(merge_id) if merge_id>1 else '',
                    rel_file)
            elif file_type=='bam' \
                and output_type=='alignments':
                input_file_param += '-ctl_filt_bam{} {} \\\n'.format(
                    bio_rep_id,                
                    rel_file)

            elif file_type=='bam' \
                and output_type=='unfiltered alignments':
                input_file_param += '-ctl_bam{} {} \\\n'.format(
                    bio_rep_id,                
                    rel_file)                
            else:
                Exception('fastq and bam input only!')

    return input_file_param.strip()

def main():
    args, ctl_exists = parse_arguments()

    mkdir_p(args.pipeline_out_root_dir)
    out_sh = os.path.join(args.pipeline_out_root_dir,
                        args.pipeline_sh_filename)
    if ctl_exists:
        map_exp_to_ctl = read_exp_to_ctl_file(args.exp_id_to_ctl_id_file)

    with open(out_sh,'w') as fp:
        sn = 0
        exp_ids = read_acc_ids_file(args.exp_acc_ids_file)
        for exp_id in exp_ids:            
            if exp_id.startswith('#'): continue
            print('==== {} ===='.format(exp_id))
            sn += 1
            exp_metadata_json_file = '{}/{}/metadata.json'.format(
                                args.exp_data_root_dir, exp_id)
            exp_metadata_org_json_file = '{}/{}/metadata.org.json'.format(
                                args.exp_data_root_dir, exp_id)
            exp_metadata_json = parse_metadata_json_file(
                exp_metadata_json_file,
                args.exp_file_type)

            species = infer_species(exp_metadata_org_json_file)

            exp_paired_end = is_paired_end(exp_metadata_org_json_file)
            if exp_paired_end:
                input_end_param = '-pe '
            else:
                input_end_param = '-se '

            if ctl_exists and exp_id in map_exp_to_ctl:
                ctl_id = map_exp_to_ctl[exp_id]
                ctl_metadata_json_file = '{}/{}/metadata.json'.format(
                                    args.ctl_data_root_dir, ctl_id)
                ctl_metadata_org_json_file = '{}/{}/metadata.org.json'.format(
                                    args.ctl_data_root_dir, ctl_id)
                ctl_metadata_json = parse_metadata_json_file(
                    ctl_metadata_json_file,
                    args.ctl_file_type if args.ctl_file_type else args.exp_file_type)
                ctl_paired_end = is_paired_end(ctl_metadata_org_json_file)
                if ctl_paired_end:
                    input_end_param += '-ctl_pe '
                else:
                    input_end_param += '-ctl_se '
            else:
                ctl_metadata_json = None
            
            input_file_param = parse_exp_metadata_json(
                exp_metadata_json, ctl_metadata_json)

            sh_item = PIPELINE_SH_ITEM_TEMPLATE.format(
                sn = sn,
                title = exp_id,
                species = species,
                input_end_param = input_end_param,
                input_file_param = input_file_param,
                pipeline_out_root_dir = args.pipeline_out_root_dir,
                pipeline_bds_script = args.pipeline_bds_script)
            fp.write('{}\n'.format(sh_item))

if __name__=='__main__':
    main()

'''
cd /oak/stanford/groups/akundaje/leepc12/run/nature_paper_recompute_xcor
python /oak/stanford/groups/akundaje/leepc12/code/ENCODE_get_ctl_from_exp/generate_pipeline_run_sh.py \
--exp-acc-ids-file /oak/stanford/groups/akundaje/leepc12/code/ENCODE_get_ctl_from_exp/exp_ids.txt \
--exp-data-root-dir /oak/stanford/groups/akundaje/leepc12/data/nature_paper_recompute_xcor/bam_exp \
--pipeline-out-root-dir . --pipeline-species hg38 --pipeline-bds-script $CODE/chipseq-pipeline/chipseq.bds \
--ctl-data-root-dir /oak/stanford/groups/akundaje/leepc12/data/nature_paper_recompute_xcor/bam_ctl \
--exp-id-to-ctl-id-file /oak/stanford/groups/akundaje/leepc12/code/ENCODE_get_ctl_from_exp/exp_ctl_ids.txt
'''