"""
test_utils.py
~~~~~~~~~~~~~
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: MIT
"""
from os.path import realpath, dirname, join

import pytest

from vep2lovd.utils import *

@pytest.fixture
def single_ped():
    return Ped(join(dirname(realpath(__file__)), "data/test.ped"))

@pytest.fixture
def trio_ped():
    return Ped(join(dirname(realpath(__file__)), "data/test_trio.ped"))

def test_numeric():
    assert make_lovd_numeric(0.3) == "0.3"
    assert make_lovd_numeric(1e-6) == "0.000001"


def test_collapse_values():
    assert collapse_values([1, 1, 1]) == "1"
    assert collapse_values([0.5, 10, 10]) == "0.5"
    assert collapse_values("a_string") == "a_string"
    assert collapse_values(1) == "1"


class TestPed(object):

    def test_validate(self, single_ped, trio_ped):
        for sample in single_ped.ped_samples:
            assert single_ped.validate(sample)
        for tr_sample in trio_ped.ped_samples:
            assert trio_ped.validate(tr_sample)

    def test_size(self, single_ped, trio_ped):
        assert len(single_ped.ped_samples) == 1
        assert len(trio_ped.ped_samples) == 3

    def test_index_sample(self, single_ped, trio_ped):
        assert len(single_ped.get_index_samples()) == 1
        assert len(trio_ped.get_index_samples()) == 1
        assert single_ped.get_index_samples()[0].individual_id == "NA12878"
        assert trio_ped.get_index_samples()[0].individual_id == "child"

    def test_is_functions(self, single_ped, trio_ped):
        assert not single_ped.is_father(single_ped.ped_samples[0])
        assert not single_ped.is_mother(single_ped.ped_samples[0])

        assert not trio_ped.is_father(trio_ped.ped_samples[0])
        assert not trio_ped.is_mother(trio_ped.ped_samples[0])

        assert trio_ped.is_father(trio_ped.ped_samples[1])
        assert not trio_ped.is_mother(trio_ped.ped_samples[1])

        assert not trio_ped.is_father(trio_ped.ped_samples[2])
        assert trio_ped.is_mother(trio_ped.ped_samples[2])

    def test_get_ped(self, single_ped, trio_ped):
        assert single_ped.get_ped("test") is None
        assert single_ped.get_ped("NA12878")

        assert trio_ped.get_ped("test") is None
        assert trio_ped.get_ped("child")
        assert trio_ped.get_ped("father")
        assert trio_ped.get_ped("mother")
