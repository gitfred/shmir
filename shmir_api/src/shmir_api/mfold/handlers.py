"""
Handlers to comunicate with mfold

"""

from flask import send_file

from shmir_api.decorators import require_json
from .mfold import mfold


@require_json(jsonify=False)
def get_mfold(data=None, **kwargs):
    filename = mfold(data)

    return send_file(filename)


get_mfold.methods = ['POST']
