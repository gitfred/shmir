#!/usr/bin/env python
import os
import logging

import json
import urllib2

from zipfile import ZipFile


def get_list(file_path, dir):
    zip_file = ZipFile(file_path, 'r')
    zip_list = zip_file.namelist()
    logging.info(zip_list)
    file_path = os.path.abspath(file_path)
    for z_file in zip_list:
        zip_file.extract(z_file, path=file_path)
    return file_path


def post_string(url, data, HEADERS):
    request = urllib2.Request(url, data, HEADERS)
    # response = urllib2.urlopen(request)
    with open("dat.zip", "w") as f:
       f.write(urllib2.urlopen(request).read())

HEADERS = {'content-type': 'application/json'}
data = json.dumps({'data': 'ACGTACGTACGTACGTACGT'})
url = 'http://127.0.0.1:5000/mfold'

post_string(url, data, HEADERS)