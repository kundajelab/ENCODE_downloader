# ENCODE downloader

This Python script downloads data files from the [ENCODE portal](https://www.encodeproject.org).

# Authentication

To download unpublished files visible to sumitters only, you need to have authentication information from the ENCODE portal.
Get your access key ID and secret key from the portal homepage menu YourID->Profile->Add Access key.

# Usage

```
usage: ENCODE downloader [-h] [--dir DIR]
                         [--file-types FILE_TYPES [FILE_TYPES ...]]
                         [--assemblies ASSEMBLIES [ASSEMBLIES ...]]
                         [--encode-access-key-id ENCODE_ACCESS_KEY_ID]
                         [--encode-secret-key ENCODE_SECRET_KEY]
                         [--ignored-accession-ids-file IGNORED_ACCESSION_IDS_FILE]
                         [--dry-run] [--dry-run-list-accession-ids]
                         [--max-download MAX_DOWNLOAD]
                         [--assembly-map ASSEMBLY_MAP [ASSEMBLY_MAP ...]]
                         [--pipeline-atac-bds-path PIPELINE_ATAC_BDS_PATH]
                         [--pipeline-chipseq-bds-path PIPELINE_CHIPSEQ_BDS_PATH]
                         [--pipeline-web-url-base PIPELINE_WEB_URL_BASE]
                         [--pipeline-encode-lab PIPELINE_ENCODE_LAB]
                         [--pipeline-encode-award PIPELINE_ENCODE_AWARD]
                         [--pipeline-encode-award-rfa PIPELINE_ENCODE_AWARD_RFA]
                         [--pipeline-encode-alias-prefix PIPELINE_ENCODE_ALIAS_PREFIX]
                         [--ignore-released | --ignore-unpublished]
                         url-or-file [url-or-file ...]

Downloads genome data files from the ENCODE portal and save them under
[WORK_DIR]/[ACCESSION_ID]/. Also generates Kundaje lab BDS pipeline shell
script for all samples (for fastq and bam only). If authentication information
(--encode-access-key-id and --encode-secret-key) is given, unpublished files
only visible to submitters with valid authentication can be downloaded.

positional arguments:
  url-or-file           List of URLs/files/accession_ids (ENCODE
                        search/experiment URL, exp. accesion ids text file or
                        exp. accession id). Make sure that URL is quotted.
                        e.g. ENCSR000ELE exp_acc_ids.txt
                        "https://www.encodeproject.org/search/?
                        type=Experiment&assay_term_name=DNase-seq &replicates.
                        library.biosample.donor.organism.scientific_name=Homo+
                        sapiens &biosample_term_name=CD14-positive+monocyte
                        &month_released=October%2C+2011".

optional arguments:
  -h, --help            show this help message and exit
  --dir DIR             [WORK_DIR] : Root directory for all downloaded genome
                        data.
  --file-types FILE_TYPES [FILE_TYPES ...]
                        List of file types to be downloaded: fastq (default),
                        bam, bed, bigWig, bigBed, ... e.g. --file-types fastq
                        bam.
  --assemblies ASSEMBLIES [ASSEMBLIES ...]
                        Assemblies (reference used for mapping) allowed. e.g.
                        --assemblies GRCh38 hg19.
  --encode-access-key-id ENCODE_ACCESS_KEY_ID
                        To download unpublished files visible to submitters
                        vith a valid authentication. Get your access key ID
                        and secret key from the portal homepage menu
                        YourID->Profile->Add Access key
  --encode-secret-key ENCODE_SECRET_KEY
                        ENCODE secret key (--encode-access-key-id must be
                        specified).
  --ignored-accession-ids-file IGNORED_ACCESSION_IDS_FILE
                        Text file with ignored accession IDs.
  --dry-run             Dry-run: downloads nothing, but generates pipeline
                        shell script.
  --dry-run-list-accession-ids
                        Dry-run: downloads nothing, but show a list of
                        accession IDs matching URL.
  --max-download MAX_DOWNLOAD
                        Maximum number of files for concurrent downloading.
                        Parallel downloading will be disabled when used with
                        --encode-access-key-id
  --assembly-map ASSEMBLY_MAP [ASSEMBLY_MAP ...]
                        List of strings to infer ENCODE assembly from species
                        name; [SPECIES_NAME]:[ASSEMBLY]. e.g. --assembly-map
                        Mus+musculus:mm10 Homo+sapiens:GRCh38
  --pipeline-atac-bds-path PIPELINE_ATAC_BDS_PATH
                        Path for atac.bds.
  --pipeline-chipseq-bds-path PIPELINE_CHIPSEQ_BDS_PATH
                        Path for chipseq.bds.
  --pipeline-web-url-base PIPELINE_WEB_URL_BASE
                        URL base for browser tracks. e.g.
                        http://mitra.stanford.edu/kundaje
  --pipeline-encode-lab PIPELINE_ENCODE_LAB
                        ENCODE lab for pipeline. e.g. /labs/anshul-kundaje/
  --pipeline-encode-award PIPELINE_ENCODE_AWARD
                        ENCODE award for pipeline. e.g. /awards/U41HG007000/
  --pipeline-encode-award-rfa PIPELINE_ENCODE_AWARD_RFA
                        ENCODE award RFA for pipeline. e.g. ENCODE3
  --pipeline-encode-alias-prefix PIPELINE_ENCODE_ALIAS_PREFIX
                        ENCODE alias prefix for pipeline (pipeline output
                        files will have aliases of [prefix]:[filename], lab
                        name is recommended). e.g. anshul-kundaje)
  --ignore-released     Ignore released data (except fastqs).
  --ignore-unpublished  Ignore unpublished data.
```

# Examples

## Using search query URL

python ENCODE_downloader.py "https://www.encodeproject.org/search/?type=Experiment&assay_title=ATAC-seq&replicates.library.biosample.life_stage=postnatal&status=released" --pipeline-encode-lab /labs/anshul-kundaje/ --pipeline-encode-award /awards/U41HG007000/ --pipeline-encode-alias-prefix anshul-kundaje --pipeline-encode-award-rfa ENCODE3 --ignore-unpublished --encode-access-key-id [ENCODE_KEY_ID] --encode-secret-key [ENCODE_PASSWORD]

## Using experiment URL

python ENCODE_downloader.py "https://www.encodeproject.org/experiments/ENCSR000ELE" --pipeline-encode-lab /labs/anshul-kundaje/ --pipeline-encode-award /awards/U41HG007000/ --pipeline-encode-alias-prefix anshul-kundaje --pipeline-encode-award-rfa ENCODE3 --ignore-unpublished --encode-access-key-id [ENCODE_KEY_ID] --encode-secret-key [ENCODE_PASSWORD]

## Using accession ids file

python ENCODE_downloader.py acc_ids.txt --pipeline-encode-lab /labs/anshul-kundaje/ --pipeline-encode-award /awards/U41HG007000/ --pipeline-encode-alias-prefix anshul-kundaje --pipeline-encode-award-rfa ENCODE3 --ignore-unpublished --encode-access-key-id [ENCODE_KEY_ID] --encode-secret-key [ENCODE_PASSWORD]

## Using accession id and ids file (mixed)

python ENCODE_downloader.py acc_ids.txt ENCSR000ELE --pipeline-encode-lab /labs/anshul-kundaje/ --pipeline-encode-award /awards/U41HG007000/ --pipeline-encode-alias-prefix anshul-kundaje --pipeline-encode-award-rfa ENCODE3 --ignore-unpublished --encode-access-key-id [ENCODE_KEY_ID] --encode-secret-key [ENCODE_PASSWORD]


