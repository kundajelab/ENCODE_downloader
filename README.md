# ENCODE downloader


It downloads experiment fastqs from the [ENCODE portal](https://www.encodeproject.org).
For the downloaded fastqs, it also creates a bash script to run genomic pipelines in Kundaje lab.

# Usage

```
ENCODE fastq downloader by Jin Lee (leepc12@gmail.com), Nov. 2016

Usage:
    python2 ENCODE_downloader.py [DIR_DOWNLOAD] [AWARD_RFA] [ASSAY_CATEGORY] [ASSAY_TITLE] [SCIENTIFIC_NAME] [REF_GENOME: optional]
Example:
    python2 ENCODE_downloader.py /scratch/data/mouse ENCODE3 DNA+accessibility ATAC-seq Mus+musculus mm10_ENCODE
    python2 ENCODE_downloader.py /scratch/data/human ENCODE3 DNA+accessibility ATAC-seq Homo+sapiens hg38_ENCODE

Do not use a space (' ') in parameters, use '+' instead.
[REF_GENOME] is optional for generating shell scripts to run BDS pipelines (using bds_scr).
```