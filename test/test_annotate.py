"""
test_annotate.py
~~~~~~~~~~~~~
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: MIT
"""

from os.path import realpath, join, dirname

import vcf
import pytest

from vep2lovd.parser import ExplodeVepVCF
from vep2lovd.utils import Ped
from vep2lovd.annotate import Annotator

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


class TestAnnotator(object):
    def test_is_refseq(self, annot, record):
        transcripts = ['ENSESTT00000020075',
                       'ENSESTT00000020087',
                       'CCDS127.1',
                       'NM_004958.3',
                       'ENSESTT00000020067',
                       'XM_005263439.1',
                       'NR_046600.1',
                       'XM_005263438.1',
                       'XM_005263440.1',
                       'XM_005263441.1']
        assert [annot.is_refseq(x) for x in record] == [False, False, False, True, False, True, False, True, True, True]

    def test_custom_annotations(self, annot, record):
        for tr in record:
            an = annot.get_custom_annotations(tr)
            assert "REF_1000G_SNP" in an
            assert "dbSNP" in an
            assert "DUPLICATED_COUNT" in an

            if an['REF_1000G_SNP'] != 'unknown':
                assert float(an['REF_1000G_SNP']) < 1.0

            if an['DUPLICATED_COUNT'] != 'unknown':
                assert int(an['DUPLICATED_COUNT']) == 3

    def test_record_values(self, annot, record):
        for tr in record:
            val = annot._get_record_values(tr)
            assert "scorePhastCons" in val

            if val['scorePhastCons'] != 'unknown':
                assert float(val['scorePhastCons']) % 1 == 0

    def test_sample_columns(self, annot, record):
        for tr in record:
            val = annot._get_sample_columns(tr)
            assert "GT_NA12878" in val
            assert "GQ_NA12878" in val
            assert "DP_NA12878" in val
            assert "DPREF_NA12878" in val
            assert "DPALT_NA12878" in val
            assert "ALTPERC_NA12878" in val
            assert "PLALT_NA12878" in val

            gt = val['GT_NA12878']
            assert gt == "0/0" or gt == "1/1" or gt == "0/1"
            gq = val['GQ_NA12878']
            assert gq == 'unknown' or (float(gq) % 1 == 0 and (0 < float(gq) < 100))
            dp = val['DP_NA12878']
            assert dp == 'unknown' or float(dp) % 1 == 0
            dpref = val['DPREF_NA12878']
            assert dpref == 'unknown' or float(dpref) % 1 == 0
            dpalt = val['DPALT_NA12878']
            assert dpalt == 'unknown' or float(dpalt) % 1 == 0
            altperc = val['ALTPERC_NA12878']
            assert altperc == 'unknown' or (0 < float(altperc) < 100)
            pltalt = val['PLALT_NA12878']
            assert pltalt == 'unknown' or float(pltalt) % 1 == 0

    def test_non_index_columns(self, annot, record):
        for tr in record:
            line = annot._get_sample_columns(tr)
            val = annot._get_sample_nonindex_columns(tr, line)
            assert val == {}

    def test_index(self, annot):
        assert annot.index == "NA12878"

    def test_distance_to_splice(self, annot, record):
        for tr in record:
            val = annot._get_distance_to_splice(tr, "DSP")
            assert "DSP" in val
            assert val["DSP"] == 'unknown' or float(val["DSP"]) % 1 == 0

    def test_bw(self, annot, record):
        bw_dir = join(dirname(realpath(__file__)), "data")
        assert annot.bigwig_directory(record[0], bw_dir) == "1"
