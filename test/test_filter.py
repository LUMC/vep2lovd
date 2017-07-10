"""
test_filter.py
~~~~~~~~~~~~~
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: MIT
"""
from os.path import realpath, join, dirname
from collections import namedtuple

import vcf
import pytest

from vep2lovd.parser import ExplodeVepVCF
from vep2lovd.utils import Ped
from vep2lovd.annotate import Annotator
from vep2lovd.filter import SingleFilter, MultiFilter

vepped = join(dirname(realpath(__file__)), "data/test.vep.vcf")
conf = join(dirname(realpath(__file__)), "data/test_annotation.json")


@pytest.fixture
def ped():
    return Ped(join(dirname(realpath(__file__)), "data/test.ped"))


@pytest.fixture
def annot(ped):
    return Annotator(conf, ped, [])


@pytest.fixture
def record():
    vep_vcf = vcf.Reader(filename=vepped)
    explod = ExplodeVepVCF(vep_vcf)
    return next(explod)


@pytest.fixture
def single_filter(annot):
    return SingleFilter(annot, ["intron_variant"], ["nc_transcript_variant"])


@pytest.fixture
def multi_filter(annot):
    return MultiFilter(annot)


class TestMultiFilter(object):
    def test_filter_xm(self, multi_filter, record):
        assert len(multi_filter.filter_xm(record)) == 1
        assert multi_filter.filter_xm(record)[0].INFO['Feature'] == "NM_004958.3"

    def test_filter_same_transcripts(self, multi_filter):
        mock_trs = ["NM_01", "NM_02", "NM_03", "NM_04", "NM_04"]
        mock_obj = namedtuple("mock_object", ["INFO"])
        mock_records = [mock_obj(INFO={"Feature": x}) for x in mock_trs]

        returned = multi_filter.filter_same_transcripts(mock_records)
        assert len(returned) == 4
        assert [x.INFO['Feature'] for x in returned] == mock_trs[:4]


class TestSingleFilter(object):
    def test_filter_consequences(self, single_filter, record):
        assert [single_filter.filter_consequences(x) for x in record] == [True, False, True, True,
                                                                          True, True, True, True, True, False]

    def test_filter_contains_consequences(self, single_filter, record):
        assert [single_filter.filter_contains_consequences(x) for x in record] == [True, True, True, True, True,
                                                                                   True, False, True, True, True]

    def test_filter_distance_to_splice(self, single_filter, record):
        assert [single_filter.filter_distance_to_splice(x, 50) for x in record] == [False, True, False, False, False,
                                                                                    False, False, False, False, True]