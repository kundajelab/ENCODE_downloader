#!/usr/bin/env python

import json
import os
import sys
import subprocess
import argparse
import signal
import collections

QUERY_URL_TEMPLATE = 'https://www.encodeproject.org/experiments/{}/?format=json'

def parse_arguments():
    parser = argparse.ArgumentParser(prog='exp_id.txt (exp_id) -> \
                        exp_to_ctl.txt (exp_id\\tctl_id)',
                        description='')
    parser.add_argument('--exp-acc-ids-file', type=str, required=True,
                            help='File with experiment accession id in each line.')
    parser.add_argument('--out-filename-exp-to-ctl', type=str, default='exp_to_ctl.txt',
                            help='exp_to_ctl.txt')
    parser.add_argument('--out-filename-ctl', type=str, default='ctl_ids.txt',
                            help='ctl_ids.txt')
    args = parser.parse_args()

    return args

def read_acc_ids(f):
    acc_ids=[]
    with open(f,'r') as fp:
        lines = fp.readlines()
        for line in lines:
            acc_ids.append(line.strip())
    return acc_ids

def get_ctl_acc_id_from_exp_acc_id(exp_acc_id):    
    json_file = '{}.json'.format(exp_acc_id)
    try:
        run_shell_cmd('wget {} -O {}'.format(
            QUERY_URL_TEMPLATE.format(exp_acc_id),
            json_file))
        json_obj = json.load(open(json_file,'r'))
        ctl_acc_ids = []
        for possible_control in json_obj["possible_controls"]:
            ctl = possible_control["@id"]
            ctl_acc_id = ctl.split('/')[2]
            ctl_acc_ids.append(ctl_acc_id)
        rm_f(json_file)
    except:
        rm_f(json_file)
        return 'NO_PERMISSION'
    return ctl_acc_ids

def rm_f(files):
    if files:
        if type(files)==list:
            run_shell_cmd('rm -f {}'.format(' '.join(files)))
        else:
            run_shell_cmd('rm -f {}'.format(files))

def run_shell_cmd(cmd): 
    try:
        p = subprocess.Popen(cmd, shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            preexec_fn=os.setsid)
        pid = p.pid
        pgid = os.getpgid(pid)
        print('run_shell_cmd: PID={}, CMD={}'.format(pid, cmd))
        ret = ''
        while True:
            line = p.stdout.readline()
            if line=='' and p.poll() is not None:
                break
            # log.debug('PID={}: {}'.format(pid,line.strip('\n')))
            print('PID={}: {}'.format(pid,line.strip('\n')))
            ret += line
        p.communicate() # wait here
        if p.returncode > 0:
            raise subprocess.CalledProcessError(
                p.returncode, cmd)
        return ret.strip('\n')
    except:
        # kill all child processes
        log.exception('Unknown exception caught. '+ \
            'Killing process group {}...'.format(pgid))
        os.killpg(pgid, signal.SIGKILL)
        p.terminate()
        raise Exception('Unknown exception caught. PID={}'.format(pid))

def main():
    args = parse_arguments()
    exp_acc_ids = read_acc_ids(args.exp_acc_ids_file)

    ctl_acc_ids = set()
    with open(args.out_filename_exp_to_ctl,'w') as fp:
        for exp_acc_id in exp_acc_ids:
            ctl_acc_id = get_ctl_acc_id_from_exp_acc_id(exp_acc_id)
            fp.write('{}\t{}\n'.format(exp_acc_id, ','.join(ctl_acc_id)))
            ctl_acc_ids.update(ctl_acc_id)

    with open(args.out_filename_ctl,'w') as fp:
        for ctl_acc_id in ctl_acc_ids:
            fp.write('{}\n'.format(ctl_acc_id))

if __name__=='__main__':
    main()
