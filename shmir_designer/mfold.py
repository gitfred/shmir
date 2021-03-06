"""
Functions for getting data from RESTful API
"""

import os
from datetime import datetime
import logging

import urllib2
import json
from zipfile import ZipFile


# URL = 'http://150.254.78.155:5000/mfold'
URL = 'http://127.0.0.1:5000/mfold'

HEADERS = {'content-type': 'application/json'}


def mfold(data=None):
    """Input: sequence string
    Output: list of files downloaded from RESTful API"""
    json_data = {'data': data}

    req = urllib2.Request(URL, json.dumps(json_data), HEADERS)

    directory = "mfold_files/"
    current_datetime = datetime.now().strftime('%H:%M:%S-%d-%m-%y')
    new_zip = os.path.join(os.path.dirname(__file__), directory,
                           'newzip_'+current_datetime+'.zip')

    if not os.path.exists(directory):
        os.makedirs(directory)

    try:
        with open(new_zip, "wb") as f:
            f.write(urllib2.urlopen(req).read())
    except urllib2.URLError:
        logging.error('Connection to mfold server refused')
        return {'error': 'Connection to mfold server refused'}
    files = get_list(new_zip, directory)
    os.remove(new_zip)
    return sorted(files)


def get_list(file_path, extract_to):
    """Function unzipps file and returns list of files"""
    with ZipFile(file_path) as zip_file:
        zip_list = zip_file.namelist()
        zip_file.extractall(path=extract_to)
        zip_list = map(lambda x: extract_to + x, zip_list)
    return zip_list
