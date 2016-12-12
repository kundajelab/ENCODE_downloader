# ENCODE downloader


It downloads experiment fastqs from the [ENCODE portal](https://www.encodeproject.org).
For the downloaded fastqs, it also creates a bash script to run genomic pipelines in Kundaje lab.

# Usage

```
usage: ENCODE fastq downloader [-h] [--encode-url-base ENCODE_URL_BASE]
                               [--ignored-accession-ids-file IGNORED_ACCESSION_IDS_FILE]
                               [--max-download MAX_DOWNLOAD]
                               [--pipeline-num-thread PIPELINE_NUM_THREAD]
                               [--pipeline-run-dir PIPELINE_RUN_DIR]
                               [--pipeline-data-dir PIPELINE_DATA_DIR]
                               [--pipeline-script PIPELINE_SCRIPT]
                               [--pipeline-web-url-base PIPELINE_WEB_URL_BASE]
                               [--pipeline-ref-genome PIPELINE_REF_GENOME]
                               [--pipeline-encode-lab PIPELINE_ENCODE_LAB]
                               [--pipeline-encode-award PIPELINE_ENCODE_AWARD]
                               [--pipeline-encode-assembly PIPELINE_ENCODE_ASSEMBLY]
                               [--pipeline-encode-alias-prefix PIPELINE_ENCODE_ALIAS_PREFIX]
                               dir-download award-rfa assay-category
                               assay-title scientific-name

Download fastqs from the ENCODE portal and generate Kundaje lab BDS pipeline
shell script for all samples.

positional arguments:
  dir-download          Root directory to save downloaded fastqs
  award-rfa             Award RFA (e.g. ENCODE3)
  assay-category        Assay category (e.g. DNA+accessibility)
  assay-title           Assay title (e.g. ATAC-seq)
  scientific-name       Scientific name for genome (e.g. Mus+musculus,
                        Homo+sapiens)

optional arguments:
  -h, --help            show this help message and exit
  --encode-url-base ENCODE_URL_BASE
                        URL base for the ENCODE portal (e.g.
                        https://www.encodeproject.org
  --ignored-accession-ids-file IGNORED_ACCESSION_IDS_FILE
                        Accession IDs in this text file will be ignored. (1
                        acc. ID per line)
  --max-download MAX_DOWNLOAD
                        Maximum number of fastqs for concurrent downloading
  --pipeline-num-thread PIPELINE_NUM_THREAD
                        Number of threads for each pipeline
  --pipeline-run-dir PIPELINE_RUN_DIR
                        Root directory of outputs for pipeline
  --pipeline-data-dir PIPELINE_DATA_DIR
                        Root directory of data for pipeline (recommended to
                        match with the first parameter [DIR_DOWNLOAD])
  --pipeline-script PIPELINE_SCRIPT
                        BDS pipeline script (e.g. /path/to/atac.bds)
  --pipeline-web-url-base PIPELINE_WEB_URL_BASE
                        URL base for browser tracks (e.g.
                        http://mitra.stanford.edu/kundaje
  --pipeline-ref-genome PIPELINE_REF_GENOME
                        Reference genome name for pipeline (e.g. hg38_ENCODE3,
                        mm10_ENCODE3, hg19, mm9)
  --pipeline-encode-lab PIPELINE_ENCODE_LAB
                        ENCODE lab for pipeline (e.g. /labs/anshul-kundaje/)
  --pipeline-encode-award PIPELINE_ENCODE_AWARD
                        ENCODE award for pipeline (e.g. /awards/U41HG007000/)
  --pipeline-encode-assembly PIPELINE_ENCODE_ASSEMBLY
                        ENCODE assembly (ref. genome name in ENCODE database)
                        for pipeline (e.g. hg38 (x), GRCh38 (o)
  --pipeline-encode-alias-prefix PIPELINE_ENCODE_ALIAS_PREFIX
                        ENCODE alias prefix for pipeline (pipeline output
                        files will have aliases of [prefix].[filename])
```

# Examples

ENCODE3 DNA+accessibility ATAC-seq Homo+sapiens
```
python2 ENCODE_downloader.py /srv/scratch/shared/surya/leepc12/data ENCODE3 DNA+accessibility ATAC-seq Homo+sapiens \
--pipeline-ref-genome hg38_ENCODE3 \
--pipeline-script '${PIPELINE_CODE}/atac_dnase_pipelines/atac.bds' \
--pipeline-num-thread 3 \
--pipeline-web-url-base http://mitra.stanford.edu/kundaje \
--pipeline-encode-lab /labs/anshul-kundaje/ \
--pipeline-encode-award /awards/U41HG007000/ \
--pipeline-encode-assembly GRCh38 \
--pipeline-encode-alias-prefix anshul-kundaje \
--ignored-accession-ids-file ignored_accession_ids.txt
```

ENCODE3 DNA+accessibility ATAC-seq Mus+musculus
```
python2 ENCODE_downloader.py /srv/scratch/shared/surya/leepc12/data ENCODE3 DNA+accessibility ATAC-seq Mus+musculus \
--pipeline-ref-genome mm10_ENCODE3 \
--pipeline-script '${PIPELINE_CODE}/atac_dnase_pipelines/atac.bds' \
--pipeline-num-thread 3 \
--pipeline-web-url-base http://mitra.stanford.edu/kundaje \
--pipeline-encode-lab /labs/anshul-kundaje/ \
--pipeline-encode-award /awards/U41HG007000/ \
--pipeline-encode-assembly mm10 \
--pipeline-encode-alias-prefix anshul-kundaje \
--ignored-accession-ids-file ignored_accession_ids.txt
```
