"""
vep2lovd.parser
~~~~~~~~~~~~~
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: MIT
"""
__author__ = 'Sander Bollen'

from collections import namedtuple, OrderedDict
import six

import vcf

from .utils import RecordWrapper, set_config, clean_transcript, Ped, collapse_values, _is_vcf_version_at_least_0_6_8
from .filter import SingleFilter, MultiFilter
from .annotate import Annotator
from .mapping import GVSFunction, DEFAULT_GATK_MAPPING

class ExplodeVepVCF(object):
    """
    This object is a parser for VEP
    It parses VEP formatted lines into standard vcf format
    It's supposed to function like an iterator
    It needs a VCF reader object
    """
    def __init__(self, reader):
        self.reader = reader
        _ = self.version_check()
        self.infos = reader.infos  # get header from reader
        self.current_line = None
        self.current_transcripts = []
        self.header_names, self.mod_infos = self.format_infos()

    def __iter__(self):
        return self

    def next(self):
        """
        Next function
        """
        self.current_line = self.reader.next()
        self.current_transcripts = self.get_transcripts()
        mod_tr = self.format_transcripts()
        return mod_tr

    def get_transcripts(self):
        """
        Get transcripts from a VEP-formatted VCF line
        Returns a list of transcripts as pyvcf records
        """
        csq = self.current_line.INFO['CSQ']
        records = []
        for transcript in csq:
            record = RecordWrapper(self.current_line)
            record.INFO['CSQ'] = transcript
            record.samples = self.current_line.samples
            records.append(record)
        return records

    def format_transcripts(self):
        """
        Format transcripts to normal vcf line
        """
        mod_tr = []
        for transcript in self.current_transcripts:
            csq_fields = transcript.INFO['CSQ'].split('|')
            assert len(csq_fields) == len(self.header_names), \
                "Length of header names and length of found fields don't match"
            for field, name in zip(csq_fields, self.header_names):
                transcript.INFO[name] = field
            mod_tr.append(transcript)
        return mod_tr

    def format_infos(self):
        """
        Format the info fields to something usable
        Returns the formatted info fields
        As well as the names of the headers
        """
        tmp_infos = self.infos
        vep_header = self.infos['CSQ'][3]
        names = vep_header.split(': ')[1].split('|')
        for name in names:
            info = namedtuple('Info', ['id', 'num', 'type', 'desc'])
            tmp = info(id=name, num=None, type='String', desc='A VEP annotation')
            tmp_infos[name] = tmp
        return names, tmp_infos

    def version_check(self):
        """
        Checks whether the VCF version is supported
        VCF < 4.0 and > 4.2 are not supported
        """
        metadata = self.reader.metadata
        # format version should not be empty
        version = metadata['fileformat']
        assert version is not None, "Version tag is empty!"

        # if version is numeric
        try:
            v = float(version)
        except:
            # assumes there's v<number> type of version
            v = version.split('v')[1]
            v = float(v)

        if v < 4.0 or v > 4.2:
            raise
        elif v >= 4.0 and v <= 4.2:
            return True

    # python3 compatibility
    def __next__(self):
        self.next()

    @property
    def all_info_names(self):
        return self.infos.keys()


class LOVDFile(object):
    def __init__(self, vcf_filename, ped_file, settings_json, filter_consequences=None,
                 filter_contains_consequences=None, filter_dsp=50, filter_xm=True, filter_refseq=True,
                 filter_identical=True, gatk_caller_field="GATKCaller", gatk_caller_mapping=DEFAULT_GATK_MAPPING, exclude_contigs_patterns=None):
        """
        Object representing an LOVD data file
        :param vcf_filename: Path to input VCF
        :param ped_file: Path to PED file
        :param settings_json: Path to settings JSON
        :param filter_consequences: List of consequences to filter out soft
        :param filter_contains_consequences: List of consequences to filter out hard
        :param filter_dsp: Max distance to splice for intronic variants (default = 50, set to 0 to disable filter)
        :param filter_xm: Filter out transcripts which are XM if NM exists
        (only together with filter_refseq. Default = True)
        :param filter_refseq: Filter out non-refseq transcripts (default = True)
        :param filter_identical: Filter out identical transcripts (default = True)
        :param gatk_caller_field: INFO field that contains caller information
        :param gatk_caller_mapping: dictionary for mapping between `gatk_caller_field` values and LOVD GATKCaller values
        :param exclude_contigs_patterns: regex patterns for contigs to exclude from filter_refseq and filter_xm
        :return:
        """

        self.VCFReader = ExplodeVepVCF(vcf.Reader(filename=vcf_filename))
        self.__ped = Ped(ped_file)
        self.index = self.__ped.get_index_samples()[0]
        self.config = set_config(settings_json)
        if exclude_contigs_patterns:
            self.exclude_contigs_re = exclude_contigs_patterns
        else:
            self.exclude_contigs_re = []

        self.annotator = Annotator(settings_json, self.__ped, self.exclude_contigs_re)

        if filter_xm and not filter_refseq:
            raise ValueError("Cannot set filter XM without filter refseq")

        self._filter_refseq = filter_refseq
        self._filter_xm = filter_xm
        self._filter_identical = filter_identical

        if filter_contains_consequences:
            self._have_contains_cons = True
            self._filter_contains_consequences = filter_contains_consequences
        else:
            self._have_contains_cons = False

        if filter_consequences:
            self._have_cons = True
            self._filter_consequences = filter_consequences
        else:
            self._have_cons = False

        if filter_dsp > 0:
            self._do_dsp = True
            self._filter_dsp = filter_dsp
        else:
            self._do_dsp = False

        self.single_filter = SingleFilter(self.annotator, filter_consequences, filter_contains_consequences)
        self.multi_filter = MultiFilter(self.annotator)

        self.mapper = GVSFunction()
        self.__current_transcripts = []
        self.__header_shown = False
        self.gatk_caller_field = gatk_caller_field
        self.gatk_caller_mapping = gatk_caller_mapping

    def apply_filters(self, trs):
        """
        Apply filters on a set of transcripts
        :param trs: transcripts
        :return: transcripts passing filters
        """

        if self._filter_identical:
            trs = self.multi_filter.filter_same_transcripts(trs)
        if self._filter_refseq and self._filter_xm:
            trs = self.multi_filter.filter_xm(trs)
        elif self._filter_refseq:
            trs = [x for x in trs if self.annotator.is_refseq(x)]

        poppables = set()

        for i, tr in enumerate(trs):
            tr = clean_transcript(tr)
            if self._have_cons:
                if not self.single_filter.filter_consequences(tr):
                    poppables.add(i)
            if self._have_contains_cons:
                if not self.single_filter.filter_contains_consequences(tr):
                    poppables.add(i)
            if self._do_dsp:
                if self.single_filter.filter_distance_to_splice(tr, self._filter_dsp):
                    poppables.add(i)
            if len(tr.INFO["Consequence"]) == 0:
                poppables.add(i)

        trs = [trs[x] for x in xrange(len(trs)) if x not in poppables]

        return trs

    def annotate(self, transcript):

        annotated = OrderedDict()

        transcript = clean_transcript(transcript)

        annotated['chromosome'] = transcript.CHROM
        annotated['position'] = transcript.POS
        annotated['REF'] = transcript.REF
        annotated['ALT'] = transcript.ALT[0]
        annotated['QUAL'] = transcript.QUAL

        if transcript.FILTER != [] and transcript.FILTER is not None:
            annotated['FILTERvcf'] = ','.join(transcript.FILTER)
        else:
            annotated['FILTERvcf'] = 'pass'

        try:
            caller_val = transcript.INFO[self.gatk_caller_field]
            if isinstance(caller_val, six.string_types):
                annotated['GATKCaller'] = self.gatk_caller_mapping[caller_val]
            else:
                annotated['GATKCaller'] = self.gatk_caller_mapping[','.join(caller_val)]
        except KeyError:
            annotated['GATKCaller'] = 'unknown'

        annotated['GVS'] = self.mapper.map(transcript.INFO["Consequence"])

        for column in self.VCFReader.all_info_names:
            if column not in transcript.INFO:
                annotated[column] = 'unknown'
                continue

            if column == 'CSQ':
                continue
            val = transcript.INFO[column]
            if val != '':
                annotated[column] = collapse_values(val)
            else:
                annotated[column] = 'unknown'

        annotated.update(self.annotator.get_custom_annotations(transcript))
        annotated.update(self.annotator._get_sample_columns(transcript))
        annotated.update(self.annotator._get_sample_nonindex_columns(transcript, annotated))
        annotated.update(self.annotator._get_record_values(transcript))
        annotated.update(self.annotator._get_distance_to_splice(transcript, 'DistanceToSplice'))

        return annotated

    def make_line(self, transcript):
        """
        Create LOVD line from a transcript
        :param transcript: VCFRecordWrapper transcript
        :return: string of LOVD line
        """

        line = []
        annotated_transcript = self.annotate(transcript)
        for _, v in annotated_transcript.iteritems():
            line.append(v)

        return "\t".join(map(str, line))

    def create_header(self, transcript):
        """
        Create LOVD header from transcript
        :param transcript: VCFRecordWrapper transcript
        :return: stirng of LOVD header line
        """
        header = []
        annotated_transcript = self.annotate(transcript)
        for k, _ in annotated_transcript.iteritems():
            header.append(k)

        return "\t".join(header)

    def next(self):
        while len(self.__current_transcripts) == 0:
            transcripts = self.VCFReader.next()
            self.__current_transcripts = self.apply_filters(transcripts)

        item = self.__current_transcripts[0]

        if not self.__header_shown:
            self.__header_shown = True
            return self.create_header(item)
        else:
            self.__current_transcripts.pop(0)
            return self.make_line(item)

    def __next__(self):
        self.next()

    def __iter__(self):
        return self
