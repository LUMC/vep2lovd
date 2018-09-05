"""
vep2lovd.utils
~~~~~~~~~~~~~
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Wibowo Arindrarto
:copyright: (c) 2016 Leiden University Medical Center
:license: MIT
"""

__author__ = 'Sander Bollen, Wibowo Arindrarto'

import json
from collections import namedtuple
import re
import vcf

config = {}


class RecordWrapper(object):
    """
    Wraps a record, such that we don't need to do deepcopy
    """
    def __init__(self, record):
        self.record = record
        self.INFO = self.copy(self.record.INFO)
        self.QUAL = self.record.QUAL
        self.REF = self.record.REF
        self.ALT = self.record.ALT
        self.CHROM = self.record.CHROM
        self.FILTER = self.record.FILTER
        self.samples = self.record.samples
        self.POS = self.record.POS
        self.is_snp = self.record.is_snp
        self.is_indel = self.record.is_indel

    def copy(self, org_dict):
        d = dict().fromkeys(org_dict)
        for k, v in org_dict.iteritems():
            try:
                d[k] = v.copy()
            except AttributeError:
                try:
                    d[k] = v[:]
                except TypeError:
                    d[k] = v
        return d


class PeekableIterator(object):

    """Class for creating iterators whose next value is viewable without
    consuming it."""

    def __init__(self, iterator):
        self._stash = []
        self._iterator = iter(iterator)

    def __iter__(self):
        return self

    def next(self):
        if self._stash:
            item = self._stash.pop(0)
        else:
            item = next(self._iterator)

        return item

    def peek(self):
        if self._stash:
            item = self._stash[0]
        else:
            try:
                item = next(self._iterator)
            except StopIteration:
                return
            self._stash = [item] + self._stash

        return item


class Ped(object):
    """
    Represents a PED file.
    See http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml for format
    """
    def __init__(self, filename):
        self.__filename = filename
        self.__handle = open(filename, "rb")
        self.PedSample = namedtuple('PedSample', [
            'family_id',
            'individual_id',
            'paternal_id',
            'maternal_id',
            'sex',
            'phenotype',
            'genotypes'])
        self.ped_samples = self._parse(self.__handle)
        self.__handle.close()
        if not all([self.validate(x) for x in self.ped_samples]):
            # should we make this a warning in stead?
            raise ValueError("Ped sample fail validation test")

    def __repr__(self):
        return "Ped(filename={0}, n_samples={1})".format(
            self.__filename,
            len(self.ped_samples)
        )

    def _parse(self, handle):
        """
        Parse handle to PED file to a list of PedSamples
        :param handle: open handle to PED file
        :return: list of PedSamples
        """
        pedsamples = []
        for line in handle:
            stripped = line.strip()
            contents = re.split('\s+', stripped)
            assert len(contents) >= 6, "Line {0} " \
                                       "is not a ped line".format(stripped)
            genotypes = []
            if len(contents) > 6:
                genotypes = contents[7:]
            pedsamples.append(self.PedSample(*contents[:6],
                                             genotypes=genotypes))
        return pedsamples

    def validate(self, ped_sample):
        """
        Validate ped sample
        :param ped_sample: PedSample instance
        :return: Boolean
        """
        alphanumeric = re.compile("^[a-zA-Z0-9_\-]*$")
        if not alphanumeric.match(ped_sample.family_id):
            return False
        if not alphanumeric.match(ped_sample.individual_id):
            return False
        if not alphanumeric.match(ped_sample.paternal_id):
            return False
        if not alphanumeric.match(ped_sample.maternal_id):
            return False
        try:
            pheno = int(ped_sample.phenotype)
        except ValueError:
            # in this case it's not numeric, and is always valid
            pass
        else:
            if pheno not in [-9, 0, 1, 2]:
                return False

        if self.is_father(ped_sample) and ped_sample.sex != "1":
            return False

        if self.is_mother(ped_sample) and ped_sample.sex != "2":
            return False

        return True

    def is_father(self, ped_sample):
        """
        Check whether ped_sample is a father
        :param ped_sample: PedSample instance
        :return: Boolean
        """
        father_ids = [x.paternal_id for x in self.ped_samples]
        if ped_sample.individual_id in father_ids:
            return True
        else:
            return False

    def is_mother(self, ped_sample):
        """
        Check whether ped_sample is a mother
        :param ped_sample: PedSample instance
        :return: Boolean
        """
        mother_ids = [x.maternal_id for x in self.ped_samples]
        if ped_sample.individual_id in mother_ids:
            return True
        else:
            return False

    def get_ped(self, id):
        """
        Get the PED at the id
        :param id: string of id
        :return: PedSample (None if not found)
        """
        ids = [x.individual_id for x in self.ped_samples]
        if id not in ids:
            return None
        return [x for x in self.ped_samples if x.individual_id == id][0]

    def get_index_samples(self):
        """
        Get affected samples
        :return: list of PedSample
        """
        return [x for x in self.ped_samples if x.phenotype == "2"]


def clean_transcript(transcript):
    """
    This function does some cleaning
    """
    # removing feature_elongation and feature_truncation
    formatted_consequence = transcript.INFO["Consequence"].replace("feature_elongation", "")
    formatted_consequence = formatted_consequence.replace("feature_truncation", "")
    if formatted_consequence.endswith("&"):
        formatted_consequence = formatted_consequence[:-1]

    transcript.INFO["Consequence"] = formatted_consequence
    return transcript


def set_config(config_file):
    """
    Parse the config file
    """
    with open(config_file, 'rb') as cf_file:
        return json.load(cf_file)


def make_lovd_numeric(number):
    """
    LOVD cannot handle number formats like 5e-05
    It instead wants everything like 0.00005

    This function formats a number to this format
    :param number: the number
    :return: string representation of number
    """
    if float(number).is_integer():
        return str(int(number))
    e_not = "{:3e}".format(number)
    if "-" not in e_not:
        return "{:.3f}".format(number)
    else:
        exp = int(float(e_not.split("-")[-1]))
        prev = e_not.split("-")[0].split(".")[-1]
        sig_num = 0
        for x in prev:
            if x != '0':
                sig_num += 1
            else:
                pass
        fmt_str = "{{:.{0}f}}".format(sig_num+exp)
        n_str = fmt_str.format(number)
    while n_str.endswith("0"):
        n_str = n_str.rstrip("0")
    return n_str



def collapse_values(val):
    """
    Return a single value for several possible inputs
    Input value can be a number, a list of number,
    a string, or a list of strings
    If it's a list, this will return the first item since
    LOVD cannot handle multiple values for a field
    :param val: the value(s)
    :return: collapsed value
    """

    if isinstance(val, list):
        first = val[0]
        # a little recursion
        return collapse_values(first)
    elif isinstance(val, str):
        try:
            return make_lovd_numeric(float(val))
        except ValueError:
            return val
    elif isinstance(val, (int, float)):
        return make_lovd_numeric(val)
    else:
        raise NotImplementedError


def _is_vcf_version_at_least_0_6_8():
    major, minor, patch = vcf.VERSION.split(".")
    if int(major) == 0 and int(minor) == 6 and int(patch) >= 8:
        return True
    elif int(major) == 0 and int(minor) > 6:
        return True
    elif int(major) > 0:
        return True
    return False
