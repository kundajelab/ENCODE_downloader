#!/usr/bin/env python

import json
import os
import sys
import argparse
import math
import collections

PIPELINE_SH_ITEM_TEMPLATE = '''#!/bin/bash
# SN={sn}
TITLE={title}
WORKDIR={pipeline_out_root_dir}/{title}; mkdir -p $WORKDIR; cd $WORKDIR
{bds} {pipeline_bds_script} -title {title} -species {species} \\
{input_file_param}
-nth {pipeline_nth_per_sample} {input_end_param} {pipeline_extra_param}
sleep 0.5
'''

def parse_arguments():
    parser = argparse.ArgumentParser(prog='Kundaje lab pipeline BDS shell script generator',
                        description='THIS PROGRAM DOES NOT SUPPORT genome hg19 and mm9!')
    parser.add_argument('--species',
                            help='-species for BDS pipeline. If not set, this will be infered from a metadata JSON.')
    parser.add_argument('--exp-acc-ids-file', type=str, required=True,
                            help='File with experiment accession id in each line.')
    parser.add_argument('--exp-data-root-dir', type=str, required=True,
                            help='Root directory where you downloaded experiment data.')
    parser.add_argument('--ctl-data-root-dir', type=str,
                            help='Root directory where you downloaded control data.')
    parser.add_argument('--exp-id-to-ctl-id-file', type=str,
                            help='File with exp_id[TAB]ctl_id in each line. Only for samples with controls.')
    parser.add_argument('--exp-file-type', type=str, required=True,
                            choices=['fastq','bam','filt_bam'],
                            help='Choose file_type of experiment replicates to run pipelines. \
                                Files with other types will be ignored exp. replicates.')
    parser.add_argument('--ctl-file-type', type=str,
                            choices=['fastq','bam','filt_bam'],
                            help='Optional. If not defined then --exp-file-type will be useds for controls. \
                                Choose file_type of control replicates to run pipelines. \
                                Files with other types will be ignored for controls.')
    parser.add_argument('--pipeline-bds-script', type=str, required=True,
                            help='Path for atac.bds or chipseq.bds.')
    parser.add_argument('--pipeline-extra-param', type=str, default='',
                            help='Extra parameter appended to the command line for BDS pipeline.\
                            Make sure that it is quoted.')
    parser.add_argument('--pipeline-out-root-dir', type=str, default='', 
                            help='Pipeline output root directory.')
    parser.add_argument('--pipeline-sh-filename-prefix', type=str, default="run_pipelines",
                            help='Prefix of path for .sh')
    parser.add_argument('--pipeline-cluster-engine', type=str, choices=['slurm','sge','local'], required=True,
                            help='Prefix of path for .sh')
    parser.add_argument('--pipeline-cluster-engine-slurm-partition', type=str,
                            help='Partition for SLURM.')
    parser.add_argument('--pipeline-cluster-engine-sge-queue', type=str,
                            help='Queue for SGE.')
    parser.add_argument('--pipeline-cluster-engine-sge-pe', type=str, default='shm',
                            help='Parallel environment for SGE.')
    parser.add_argument('--pipeline-nth-per-sample', type=int, default=3,
                            help='Number of threads per sample.')
    parser.add_argument('--pipeline-mem-per-sample', type=int, default=40,
                            help='Memory limit in GB per sample.')
    parser.add_argument('--pipeline-walltime-per-sample', type=int, default=47,
                            help='Walltime in hours per sample.')
    parser.add_argument('--pipeline-number-of-samples-per-sh', type=int, default=50,
                            help='Number of samples per .sh.')
    args = parser.parse_args()

    if args.ctl_data_root_dir and not args.exp_id_to_ctl_id_file or \
        not args.ctl_data_root_dir and args.exp_id_to_ctl_id_file:
        raise Exception('--ctl-data-root-dir and --exp-id-to-ctl-id-file must be defined together.')
    if args.pipeline_nth_per_sample<1:
        raise Exception('--pipeline-nth-per-sample must be >0.')
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
    assembly = json_obj['assembly']
    if 'GRCh38' in assembly: return 'hg38'
    if 'hg19' in assembly: return 'hg19'
    if 'GRCm38' in assembly or 'mm10' in assembly: return 'mm10'
    if 'mm9' in assembly: return 'mm9'
    if deep_search(json_obj, 'Homo sapiens'):
        return 'hg38'
    if deep_search(json_obj, 'Mus Musculus'):
        return 'mm10'
    raise Exception('could not find/infer species from {}'.format(
        metadata_org_json_file))

def get_contributing_file_acc_ids(metadata_org_json_file):
    with open(metadata_org_json_file,'r') as fp:
        json_obj = json.load(fp)
    result = []
    # convert /files/[file_acc_id]/ to [file_acc_id]
    for s in json_obj['contributing_files']:
        result.append(s.split('/files/')[1].strip('/'))
    # print(result)
    return result

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

        # print(file_type, output_type, bio_rep_id, pair, paired_with, rel_file)

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
            result.append((paired_with, file_type2, output_type2,
                new_bio_rep_id2, pair2, merge_id, paired_with2, rel_file2))
        else: # fastq se or other file_type
            merge_id = 1

        new_bio_rep_id = map_bio_rep_id_to_serial_bio_rep_id[bio_rep_id]
        result.append((file_acc_id, file_type, output_type, 
            new_bio_rep_id, pair, merge_id, paired_with, rel_file))
    return result

def parse_exp_metadata_json(exp, ctls, contributing_file_acc_ids):
    input_file_param = ''

    print("exps:")
    for (file_acc_id, file_type, output_type, bio_rep_id, pair, merge_id,
        paired_with, rel_file) in exp:

        print(file_acc_id, file_type, output_type, bio_rep_id, pair, merge_id,
            paired_with, rel_file)
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

    if ctls: # if there is control
        print("ctls:")
        cnt=0
        for ctl in ctls:
            for (file_acc_id, file_type, output_type, bio_rep_id, pair, merge_id,
                paired_with, rel_file) in ctl:
                if not file_acc_id in contributing_file_acc_ids:
                    continue
                cnt+=1
        if cnt==0: 
            # if corresponding control file_acc_id not found
            # then force to use 1st control set 
            skip_checking_file_acc_id = True
        else:
            skip_checking_file_acc_id = False

        for ctl in ctls:
            for (file_acc_id, file_type, output_type, bio_rep_id, pair, merge_id,
                paired_with, rel_file) in ctl:
                if not skip_checking_file_acc_id and not file_acc_id in contributing_file_acc_ids:
                    continue
                print(file_acc_id, file_type, output_type, bio_rep_id, pair, merge_id,
                    paired_with, rel_file)

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
            if skip_checking_file_acc_id:
                break
    return input_file_param.strip()

def main():
    args, ctl_exists = parse_arguments()

    mkdir_p(args.pipeline_out_root_dir)

    if ctl_exists:
        map_exp_to_ctl = read_exp_to_ctl_file(args.exp_id_to_ctl_id_file)

    sh_items = []

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

        if args.species:
            species = args.species
        else:
            species = infer_species(exp_metadata_org_json_file)

        exp_paired_end = is_paired_end(exp_metadata_org_json_file)
        if exp_paired_end:
            input_end_param = '-pe '
        else:
            input_end_param = '-se '

        ctl_metadata_jsons = []
        if ctl_exists and exp_id in map_exp_to_ctl:
            for ctl_id in map_exp_to_ctl[exp_id].split(','):
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
                contributing_file_acc_ids = get_contributing_file_acc_ids(exp_metadata_org_json_file)
                ctl_metadata_jsons.append(ctl_metadata_json)
        else:
            ctl_metadata_json = None
            contributing_file_acc_ids = []
        
        input_file_param = parse_exp_metadata_json(
            exp_metadata_json, ctl_metadata_jsons, contributing_file_acc_ids)

        sh_item = PIPELINE_SH_ITEM_TEMPLATE.format(
            sn = sn,            
            title = exp_id,
            bds = 'bds_scr {}'.format(exp_id) if args.pipeline_cluster_engine=='local' else 'bds',
            species = species,
            input_end_param = input_end_param,
            input_file_param = input_file_param,
            pipeline_out_root_dir = args.pipeline_out_root_dir,
            pipeline_bds_script = args.pipeline_bds_script,
            pipeline_nth_per_sample = args.pipeline_nth_per_sample,
            pipeline_extra_param = '-system local ' + args.pipeline_extra_param)
        sh_items.append((exp_id, sh_item))
    
    master_sh_prefix = os.path.join(args.pipeline_out_root_dir,
                    args.pipeline_sh_filename_prefix)
    for i in range(int(math.ceil(len(sh_items)/float(args.pipeline_number_of_samples_per_sh)))):
        start = i*args.pipeline_number_of_samples_per_sh
        end = min(len(sh_items), (i+1)*args.pipeline_number_of_samples_per_sh)
        lines_in_master_sh = ''
        # write sh for individual sample
        for j in range(start,end):
            exp_id, sh_item = sh_items[j]
            sample_sh_prefix = os.path.join(args.pipeline_out_root_dir, format(exp_id))
            sample_sh = '{}.sh'.format(sample_sh_prefix)
            with open(sample_sh,'w') as fp:
                fp.write(sh_item)
            sample_out_dir = os.path.join(args.pipeline_out_root_dir, exp_id)
            mkdir_p(sample_out_dir)
            o = os.path.join(sample_out_dir, 'out.log')
            e = o
            lines_in_master_sh += 'echo "SN={} EXP_ID={}"\n'.format(j+1,exp_id)
            if args.pipeline_cluster_engine=='slurm':
                line = 'sbatch -J {} -o {} -e {} --export=ALL -n 1 --ntasks-per-node=1 --cpus-per-task={} '
                line += '--mem {}G -t {} -p {} {}'
                line = line.format(
                    exp_id,
                    o,
                    e,                    
                    args.pipeline_nth_per_sample,
                    int(args.pipeline_mem_per_sample), # GB
                    args.pipeline_walltime_per_sample*60, # hours
                    args.pipeline_cluster_engine_slurm_partition,
                    sample_sh)
            elif args.pipeline_cluster_engine=='sge':
                line = 'qsub -o {} -e {} -V -pe {} {} '
                line += '-l h_vmem={}G,s_vmem={}G,h_rt={}:00:00,s_rt={}:00:00 -q {} {}'
                line = line.format(
                    o,
                    e,
                    args.pipeline_cluster_engine_sge_pe,
                    args.pipeline_nth_per_sample,
                    args.pipeline_mem_per_sample,
                    args.pipeline_mem_per_sample,
                    args.pipeline_walltime_per_sample,
                    args.pipeline_walltime_per_sample,
                    args.pipeline_cluster_engine_sge_queue,
                    sample_sh)
            else:
                line = 'bash {}'.format(sample_sh)
            lines_in_master_sh += '{}\nsleep 5\n\n'.format(line)

        # write master runner sh for group of sample .sh
        with open('{prefix}.{start:04d}-{end:04d}.sh'.format(
                    prefix = master_sh_prefix,
                    start = start+1,
                    end = end),'w') as fp:
            fp.write(lines_in_master_sh)
            # fp.write('\n'.join(sh_items[start:end]))
    
if __name__=='__main__':
    main()
