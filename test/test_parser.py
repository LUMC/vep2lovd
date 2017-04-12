"""
test_parser.py
~~~~~~~~~~~~~
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: MIT
"""

from os.path import realpath, join, dirname

import vcf
import pytest

from vep2lovd.parser import ExplodeVepVCF, LOVDFile

unvepped = join(dirname(realpath(__file__)), "data/test.vcf")
vepped = join(dirname(realpath(__file__)), "data/test.vep.vcf")
conf = join(dirname(realpath(__file__)), "data/test_annotation.json")
ped = join(dirname(realpath(__file__)), "data/test.ped")


class TestExplodeVepVcf(object):

    def test_unvepped_exception(self):
        un_reader = vcf.Reader(filename=unvepped)
        with pytest.raises(KeyError):
            _ = ExplodeVepVCF(un_reader)

    def test_header(self):
        vep_reader = vcf.Reader(filename=vepped)
        explod = ExplodeVepVCF(vep_reader)
        assert len(explod.header_names) == 30

    def test_has_all_fields(self):
        vep_reader = vcf.Reader(filename=vepped)
        explod = ExplodeVepVCF(vep_reader)
        for explosion in explod:
            for record in explosion:
                for h in explod.header_names:
                    assert h in record.INFO

    def test_info_fields(self):
        vep_reader = vcf.Reader(filename=vepped)
        explod = ExplodeVepVCF(vep_reader)
        expected_infos = [
            "AC", "AF", "AN", "BaseQRankSum", "DB", "DP", "DS", "Dels", "END", "FS", "HaplotypeScore",
            "InbreedingCoeff", "MLEAC", "MLEAF", "MQ", "MQ0", "MQRankSum", "NEGATIVE_TRAIN_SITE", "POSITIVE_TRAIN_SITE",
            "QD", "RPA", "RU", "ReadPosRankSum", "STR", "VQSLOD", "culprit", "ClippingRankSum", "GATKCaller", "Allele",
            "Gene", "Feature", "Feature_type", "Consequence", "cDNA_position", "CDS_position", "Protein_position",
            "Amino_acids", "Codons", "Existing_variation", "DISTANCE", "STRAND", "HGVSc", "HGVSp", "PUBMED", "PolyPhen",
            "ALLELE_NUM", "CLIN_SIG", "EXON", "INTRON", "AA_MAF", "EA_MAF", "GMAF", "SYMBOL", "SYMBOL_SOURCE",
            "AFR_MAF", "AMR_MAF", "ASN_MAF", "EUR_MAF", "CSQ", "scorePhastCons"
        ]
        assert sorted(expected_infos) == sorted(explod.all_info_names)


class TestLOVDFile(object):

    def test_conversion(self):
        lovd = LOVDFile(vepped, ped, conf, filter_refseq=False, filter_xm=False)
        i = 0
        for _ in enumerate(lovd):
            i += 1
        assert i == 284

        with_filter = LOVDFile(vepped, ped, conf)
        o = 0
        for _ in enumerate(with_filter):
            o += 1
        assert o == 51


