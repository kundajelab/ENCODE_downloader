#!/usr/bin/env python2

# import encodedcc
import requests
import urllib2
import json
import os

# key = encodedcc.ENC_Key("/users/leepc12/ENCODE_keypairs/keypairs.json", "production")
# connection = encodedcc.ENC_Connection(key)
HEADERS = {'accept': 'application/json'}
# encode_search_url = 'https://www.encodeproject.org/search/?type=Experiment&assay_title=ATAC-seq&award.project=ENCODE&lab.title=Barbara+Wold%2C+Caltech'
# search_data = requests.get(encode_search_url, headers=HEADERS, auth=connection.auth)

encode_search_url = "https://www.encodeproject.org/experiments/ENCSR229QKB"
# r = requests.get(encode_search_url, headers=HEADERS, auth=("VTCUINGZ","bem5gvvrgdqo3j37"))
r = requests.get(encode_search_url, auth=("VTCUINGZ","bem5gvvrgdqo3j37"))
print r.text
# print(search_data.text)

# r = requests.post(encode_search_url, headers=HEADERS, data="VTCUINGZ:bem5gvvrgdqo3j37")
# print r.text

# r = urllib2.urlopen(encode_search_url, auth=("VTCUINGZ","bem5gvvrgdqo3j37"))
# data = json.load(r)
# print(data)

# def DownloadFile(url):
#     local_filename = url.split('/')[-1]
#     # r = requests.get(url, headers=HEADERS, auth=connection.auth)
#     r = requests.get(url, auth=("VTCUINGZ","bem5gvvrgdqo3j37"))
#     with open(local_filename, 'wb') as f:
#         for chunk in r.iter_content(chunk_size=1024): 
#             if chunk: # filter out keep-alive new chunks
#                 f.write(chunk)
#     return 

# # DownloadFile("https://www.encodeproject.org/files/ENCFF217JOA/@@download/ENCFF217JOA.fastq.gz")
# DownloadFile("https://www.encodeproject.org/files/ENCFF196YIG/@@download/ENCFF196YIG.bigWig")


