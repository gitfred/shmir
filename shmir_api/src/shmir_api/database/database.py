"""
Module for communication with database
"""
import re
from flask import g

from sqlsoup import SQLSoup
from sqlalchemy import func
from shmir_api.settings import DB_NAME, DB_USER, DB_PASS, DB_HOST, DB_PORT


class SerializeMixin:
    """
    Mixin to serialize attributes which do not start with _ or are not id
    """
    def serialize(self):
        return {key: value for key, value in vars(self).items()
                if not key.startswith("_") and key != "id"}


def get_db():
    """
    Global connector variable
    """
    db = getattr(g, '_database', None)
    if db is None:
        conn_str = "postgresql+psycopg2://{user}:{password}@{host}:{port}/{dbname}"
        fconn = conn_str.format(dbname=DB_NAME, user=DB_USER,
                                password=DB_PASS, host=DB_HOST, port=DB_PORT)
        db = g._database = SQLSoup(fconn, SerializeMixin)

    return db


def disconnect():
    """
    Global disconnector
    """
    db = getattr(g, '_database', None)
    if db is not None:
        db.connection().close()


def serialized_all_by_query(query):
    """
    Function which returns serialized query
    """
    return [elem.serialize() for elem in query.all()]


def backbone_get_all():
    """
    Function which gets all serialized Backbones in database
    """
    db = get_db()
    return serialized_all_by_query(db.backbone)


def backbone_get_by_name(name):
    """
    Function which gets one serialized Backbone by name
    """
    db = get_db()
    data = db.backbone.filter(func.lower(db.backbone.name) ==
                              func.lower(name)).first()

    return data.serialize() if data else {}


def backbone_get_by_miRNA_s(letters):
    """
    Function which gets serialized Backbones having first two letters of
    miRNA_s same as letters given (first two nucleotides of siRNA strand)
    """
    db = get_db()
    query = db.backbone.filter(db.backbone.miRNA_s.
                               like("{}%".format(letters.upper())))

    return serialized_all_by_query(query)


def immuno_get_all():
    """
    Function which gets all serialized immuno sequences in database
    """
    db = get_db()
    return serialized_all_by_query(db.immuno)


def create_regular(seq):
    """Function for generating regular expresions for given miRNA sequence
    according to the schema below:

    example miRNA sequence: UGUAAACAUCCUCGACUGGAAG
    U... (weight 1): the first nucleotide
    U...G (weight 2): the first and the last nucleotides
    UG... (weight 2): two first nucleotides
    UG...G (weight 3): two first and the last nucleotides
    UG...AG (weight 4): two first and two last nucleotides
    """
    seq = seq.upper()
    acids = '[UTGCA]' # order is important: U should always be next to T
    generic = r'{1}{2}[UTGCA]{{{0}}}{3}{4}'
    begin = [
        {
            'acid': re.sub('[UT]', '[UT]', letter),
            'excluded': acids.replace('UT' if letter in 'UT' else letter, '')
        }
        for letter in seq[:2]]

    end = [
        {
            'acid': re.sub('[UT]', '[UT]', letter),
            'excluded': acids.replace('UT' if letter in 'UT' else letter, '')
        }
        for letter in seq[-2:]]

    regexp1 = [generic.format(i, begin[0]['acid'], begin[1]['excluded'],
                              end[0]['excluded'], end[1]['excluded'])
               for i in range(15, 18)]

    regexp2 = [generic.format(i, begin[0]['acid'], begin[1]['excluded'],
                              end[0]['excluded'], end[1]['acid'])
               for i in range(15, 18)]

    regexp2_ = [generic.format(i, begin[0]['acid'], begin[1]['acid'],
                               end[0]['excluded'], end[1]['excluded'])
                for i in range(15, 18)]

    regexp3 = [generic.format(i, begin[0]['acid'], begin[1]['acid'],
                              end[0]['excluded'], end[1]['acid'])
               for i in range(15, 18)]

    regexp4 = [generic.format(i, begin[0]['acid'], begin[1]['acid'],
                              end[0]['acid'], end[1]['acid'])
               for i in range(15, 18)]

    return {
        1: regexp1,
        2: regexp2 + regexp2_,
        3: regexp3,
        4: regexp4,
        }
