# ENCODE downloader


It downloads experiment fastqs from the [ENCODE portal](https://www.encodeproject.org).
For the downloaded fastqs, it also creates a bash script to run genomic pipelines in Kundaje lab.

# Usage

```
usage: ENCODE fastq downloader [-h]
                               [--ignored-accession-ids-file IGNORED_ACCESSION_IDS_FILE]
                               [--max-download MAX_DOWNLOAD]
                               [--ref-genome REF_GENOME]
                               [--num-thread-pipeline NUM_THREAD_PIPELINE]
                               [--dir-pipeline-run DIR_PIPELINE_RUN]
                               [--dir-pipeline-data DIR_PIPELINE_DATA]
                               [--bds-pipeline-script BDS_PIPELINE_SCRIPT]
                               [--web-url-base WEB_URL_BASE]
                               [--encode-url-base ENCODE_URL_BASE]
                               dir-download award-rfa assay-category
                               assay-title scientific-name

Download fastqs from the ENCODE portal.

positional arguments:
  dir-download          Root directory to save downloaded fastqs
  award-rfa             Award RFA (e.g. ENCODE3)
  assay-category        Assay category (e.g. DNA+accessibility)
  assay-title           Assay title (e.g. ATAC-seq)
  scientific-name       Scientific name for genome (e.g. Mus+musculus,
                        Homo+sapiens)

optional arguments:
  -h, --help            show this help message and exit
  --ignored-accession-ids-file IGNORED_ACCESSION_IDS_FILE
                        Accession IDs in this text file will be ignored. (1
                        acc. ID per line)
  --max-download MAX_DOWNLOAD
                        Maximum number of fastqs for concurrent downloading
  --ref-genome REF_GENOME
                        Reference genome name for pipeline (e.g. hg38_ENCODE3,
                        mm10_ENCODE3, hg19, mm9)
  --num-thread-pipeline NUM_THREAD_PIPELINE
                        Number of thread for each pipeline
  --dir-pipeline-run DIR_PIPELINE_RUN
                        Root directory of outputs for pipeline
  --dir-pipeline-data DIR_PIPELINE_DATA
                        Root directory of data for pipeline
  --bds-pipeline-script BDS_PIPELINE_SCRIPT
                        BDS pipeline script (e.g. /path/to/atac.bds)
  --web-url-base WEB_URL_BASE
                        URL base for browser tracks (e.g.
                        http://mitra.stanford.edu/kundaje
  --encode-url-base ENCODE_URL_BASE
                        URL base for the ENCODE portal (e.g.
                        https://www.encodeproject.org
```

# Examples

ENCODE3 DNA+accessibility ATAC-seq Homo+sapiens
```
python2 ENCODE_downloader.py /srv/scratch/shared/surya/leepc12/data ENCODE3 DNA+accessibility ATAC-seq Homo+sapiens \
--ref-genome hg38_ENCODE3 --bds-pipeline-script '${DIR_PIPELINE_CODE}/atac_dnase_pipelines/atac.bds' --num-thread-pipeline 3 \
--web-url-base http://mitra.stanford.edu/kundaje --ignored-accession-ids-file ignored_accession_ids.txt
```

ENCODE3 DNA+accessibility ATAC-seq Mus+musculus
```
python2 ENCODE_downloader.py /srv/scratch/shared/surya/leepc12/data ENCODE3 DNA+accessibility ATAC-seq Mus+musculus \
--ref-genome mm10_ENCODE3 --bds-pipeline-script '${DIR_PIPELINE_CODE}/atac_dnase_pipelines/atac.bds' --num-thread-pipeline 3 \
--web-url-base http://mitra.stanford.edu/kundaje --ignored-accession-ids-file ignored_accession_ids.txt
```
