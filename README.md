# ENCODE downloader

This Python script downloads data files from the [ENCODE portal](https://www.encodeproject.org).

# Supported data types

You can download any data type of data from the portal. For example to download all FASTQs and unfiltered BAMs:
```
$ python encode_downloader.py --file-types fastq "bam:unfiltered alignments" ...
```
For example to download all p-value signal bigwigs:
```
$ python encode_downloader.py --file-types "bigWig:signal p-value" ...
```

# Authentication

To download unpublished files visible to sumitters only, you need to have authentication information from the ENCODE portal.
Get your access key ID and secret key from the portal homepage menu YourID->Profile->Add Access key.

```
$ python encode_downloader.py --encode-access-key-id [ENCODE_ACCESS_KEY_ID] --encode-secret-key [ENCODE_SECRET_KEY] ...
```

# Usage

```
$ python encode_downloader.py -h
```

# Generating BDS pipeline script

After you download data files you need to process them with pipelines. `generate_pipeline_run_sh.py` generates a shell script `run_pipelines.sh` to run Kundaje lab's BDS pipelines.

* [ATAC-Seq/DNase-Seq pipeline](https://github.com/kundajelab/atac_dnase_pipelines)
* [TF/Histone ChIP-Seq pipeline](https://github.com/kundajelab/chipseq_pipeline)

`[EXP_ACC_IDS_TXT]` is a text file with `[experiment_accession_id]` in each line. Specify root directory of experiment data files `[EXP_DATA_ROOT_DIR]` that you downloaded from ENCODE portal by `encode_downloader.py`.

To get a shell script to run pipelines without controls:
```
$ python generate_pipeline_run_sh.py --exp-acc-ids-file [EXP_ACC_IDS_TXT] --exp-data-root-dir [EXP_DATA_ROOT_DIR] --pipeline-bds-script [BDS_FILE_PATH; {chipsqe.bds, atac.bds}] --file-type-to-run-pipeline [FILE_TYPE; {fastq,bam,filt_bam}]
```

To get a shell script To run pipelines with controls. `[CTL_ACC_IDS_TXT]` is a TSV file with `[experiment_accession_id]\t[control_accession_id]`. Specify root directory of control data files `[CTL_DATA_ROOT_DIR]` that you downloaded from ENCODE portal by `encode_downloader.py`.
```
$ python generate_pipeline_run_sh.py --exp-acc-ids-file [EXP_ACC_IDS_TXT] --exp-data-root-dir [EXP_DATA_ROOT_DIR] --exp-id-to-ctl-id-file exp_to_ctl.txt --ctl-data-root-dir [CTL_DATA_ROOT_DIR] --pipeline-bds-script [BDS_FILE_PATH; chipsqe.bds or atac.bds] --file-type-to-run-pipeline [FILE_TYPE; {fastq,bam,filt_bam}]
```

# Requirements

* Python requests
```
$ pip install requests
```

# Examples

### Using search query URL

```
$ python encode_downloader.py "https://www.encodeproject.org/search/?type=Experiment&assay_title=ATAC-seq&replicates.library.biosample.life_stage=postnatal&status=released" --ignore-unpublished --file-types fastq "bam:unfiltered alignments"
```

### Using experiment URL

```
$ python encode_downloader.py "https://www.encodeproject.org/experiments/ENCSR000ELE" --ignore-unpublished --encode-access-key-id [ENCODE_KEY_ID] --encode-secret-key [ENCODE_PASSWORD]
```

### Using accession ids file

```
$ python encode_downloader.py acc_ids.txt --ignore-unpublished --encode-access-key-id [ENCODE_KEY_ID] --encode-secret-key [ENCODE_PASSWORD] --file-types "bigWig:signal p-value"
```

### Using accession id and ids file (mixed)

```
$ python encode_downloader.py acc_ids.txt ENCSR000ELE --ignore-unpublished --encode-access-key-id [ENCODE_KEY_ID] --encode-secret-key [ENCODE_PASSWORD]
```