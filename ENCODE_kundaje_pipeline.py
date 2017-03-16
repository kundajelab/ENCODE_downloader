import collections

class PipelineShellScript(object):    
    def __init__(self, lab, alias_prefix, award,
            data_root_dir, web_url_base, chipseq_bds_path, atac_bds_path):
        self.lab = lab
        self.alias_prefix = alias_prefix
        self.award = award
        self.web_url_base = web_url_base
        self.data_root_dir = data_root_dir
        self.samples = []
        self.chipseq_bds_path = chipseq_bds_path
        self.atac_bds_path = atac_bds_path

    # files = dict(key:file_type(fastq,bam,...), val:(bio_rep_id,tech_rep_id,pair_id,))
    def add_sample(self, title, assembly, assay_title, assay_category, award_rfa, files):
        self.samples.append( PipelineSample(self, title, assembly, assay_title, assay_category, award_rfa, files) )

    def __str__(self):
        contents = '#!/bin/bash\n\nPIPELINE_RUN_DIR=$PWD\n'.format(self.data_root_dir)
        contents += 'PIPELINE_DATA_DIR={}\n'.format(self.data_root_dir)
        if self.web_url_base:
            contents += 'PIPELINE_WEB_URL_BASE={}\n'.format(self.web_url_base)
        contents += '\n'
        for sample in self.samples:
            contents += str(sample)
        return contents

    def write_to_file(self, fname='run_pipelines.sh'):
        with open(fname,'w') as f:
            f.write(str(self))

class PipelineSample(object):
    counter = 0
    # acc_id = ENCODE accession id
    def __init__(self, sh, acc_id, assembly, assay_title, assay_category, award_rfa, files):
        PipelineSample.counter += 1
        self.sn = PipelineSample.counter
        self.sh = sh
        self.acc_id = acc_id
        self.assembly = assembly
        self.assay_title = assay_title
        self.assay_category = assay_category
        self.award_rfa = award_rfa
        self.species = self.__get_species_from_assembly()
        self.script = self.__get_script_from_assay_title()
        self.extra_param = self.__get_extra_param()
        self.files = files

    def __get_species_from_assembly(self):
        if self.assembly.startswith('GRCh38') or self.assembly.startswith('hg38'):
            return 'hg38_ENCODE'
        elif self.assembly.startswith('hg19'):
            return 'hg19'
        elif self.assembly.startswith('mm10'):
            return 'mm10_ENCODE'
        elif self.assembly.startswith('mm9'):
            return 'mm9'
        return ''

    def __get_script_from_assay_title(self):
        if self.assay_title in ['DNase-seq','ATAC-seq']:
            return self.sh.atac_bds_path
        elif self.assay_title in ['ChIP-seq']:
            return self.sh.chipseq_bds_path
        return ''

    def __get_extra_param(self):
        ret = ''
        if self.assay_title == 'DNase-seq':
            ret += '-dnase_seq '
        if self.award_rfa == 'ENCODE3':
            ret += '-ENCODE3 '
        return ret

    def __get_pipeline_script(self):
        if self.assay_title =='DNase-seq':
            return self.sh.atac_bds_path + ' -dnase_seq'
        elif self.assay_title =='ATAC-seq':
            return self.sh.atac_bds_path
        elif self.assay_title =='ChIP-seq':
            return self.sh.chipseq_bds_path
        return ret

    def __has_file_type(self,file_type):
        return sum([self.files[file_acc_id][0]==file_type for file_acc_id in self.files])

    def __get_file_types(self):
        return list(set([self.files[file_acc_id][0] for file_acc_id in self.files]))

    def __str__(self):
        # header
        ret = '# SN={}, EXP. ACCESSION ID={}\n'.format(self.sn,self.acc_id)
        if not self.species:
            ret += '# Unsupported species: {}\n'.format(self.assembly)
        elif not self.script:
            ret += '# Unsupported assay title: {}\n'.format(self.assay_title)
        elif not self.__has_file_type('fastq') and not self.__has_file_type('bam'):
            ret += '# Unsupported file type: {}\n'.format(self.__get_file_types())
        else:
            # when multiple type of files exist (fastq and bam), fastq has priority
            if self.__has_file_type('bam'):
                reserved_file_type = 'bam'
            elif self.__has_file_type('fastq'):
                reserved_file_type = 'fastq'
            # mkdir
            ret += 'WORKDIR=${PIPELINE_RUN_DIR}/%s; mkdir -p $WORKDIR; cd $WORKDIR\n' % (self.acc_id,)
            # files (key: file_accession_id, val:file_type, ...)
            cnt_fastq_to_be_pooled = collections.defaultdict(int)
            already_done = []
            param_file = ''
            for file_acc_id in self.files:
                file_type, status, bio_rep_id, pair, paired_with, rel_file = self.files[file_acc_id]
                if rel_file in already_done:
                    continue                
                cnt_fastq_to_be_pooled[bio_rep_id] += 1
                # suffix for files to be pooled
                param_endedness = ('-pe' if paired_with else '-se' ) + str(bio_rep_id)
                if cnt_fastq_to_be_pooled[bio_rep_id] > 1:
                    suffix = '_'+str(cnt_fastq_to_be_pooled[bio_rep_id])
                    suffix2 = ':'+str(cnt_fastq_to_be_pooled[bio_rep_id])
                    param_endedness = ''
                else:
                    suffix = ''
                    suffix2 = ''
                if file_type==reserved_file_type:
                    if file_type == 'fastq':
                        if paired_with:
                            _, _, _, pair2, _, rel_file2 = self.files[paired_with]
                            ret += 'FASTQ{}_{}{}={}\n'.format(bio_rep_id, pair, suffix, rel_file)
                            ret += 'FASTQ{}_{}{}={}\n'.format(bio_rep_id, pair2, suffix, rel_file2)
                            param_file += ' {} -fastq{}_{}{} $FASTQ{}_{}{} -fastq{}_{}{} $FASTQ{}_{}{}' \
                                               .format( param_endedness, bio_rep_id, pair, suffix2, bio_rep_id, pair, suffix, \
                                                               bio_rep_id, pair2, suffix2, bio_rep_id, pair2, suffix )
                            already_done.append(rel_file2)
                        else:
                            ret += 'FASTQ{}{}={}\n'.format(bio_rep_id, suffix, rel_file)
                            param_file += ' {} -fastq{}{} $FASTQ{}{}' \
                                               .format(param_endedness, bio_rep_id, suffix2, bio_rep_id, suffix)
                    elif file_type=='bam':
                        if status != 'released': 
                            print('Warning: ignored unpublished bam (file_accession_id: {}, status: {}) in pipeline shell script'.format(
                                file_acc_id,status))
                            continue
                        ret += 'BAM{}={}\n'.format(bio_rep_id, rel_file)
                        param_file += ' {} -filt_bam{} $BAM{}' \
                                           .format(param_endedness, bio_rep_id, bio_rep_id)
                already_done.append(rel_file)

            # BDS command line and parameters
            param = ('-title {} -species {} {} '+\
                    '-ENCODE_accession {} '+\
                    '-ENCODE_assay_title {} '+\
                    '-ENCODE_assembly {} ').format(
                        self.acc_id, self.species, param_file,
                        self.acc_id,
                        self.assay_title,
                        self.assembly)    

            if self.award_rfa:
                param += '-ENCODE_award_rfa {} '.format(self.award_rfa)
            if self.assay_category:
                param += '-ENCODE_assay_category {} '.format(self.assay_category)
            if self.sh.web_url_base:
                param += '-url_base {} '.format(self.sh.web_url_base+'/'+self.acc_id+'/out')
            if self.sh.lab:
                param += '-ENCODE_lab {} '.format(self.sh.lab)
            if self.sh.award: 
                param += '-ENCODE_award {} '.format(self.sh.award)
            if self.sh.alias_prefix: 
                param += '-ENCODE_alias_prefix {} '.format(self.sh.alias_prefix)
            param += self.extra_param
            pipeline_script = self.__get_pipeline_script()

            ret += 'bds_scr $TITLE {} {}\nsleep 0.5\n\n'.format(pipeline_script, param)
        # ret += '\n\n'
        return ret
