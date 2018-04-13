"""
vep2lovd.annotate
~~~~~~~~~~~~~
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Wibowo Arindrarto
:copyright: (c) 2016 Leiden University Medical Center
:license: MIT
"""

__author__ = 'Sander Bollen, Wibowo Arindrarto'

import os
import sys

import vcf
import pyBigWig
from pysam import Tabixfile
import hgvs.parser
from hgvs.exceptions import HGVSParseError

from .utils import PeekableIterator, set_config, \
    collapse_values, _is_vcf_version_at_least_0_6_8

COLS_RECORD = [
    'scorePhastCons'
]

COLS_PER_SAMPLE = ['GT', 'GQ', 'DP', 'DPREF', 'DPALT',
                   'ALTPERC', 'PLALT']

COLS_PER_SAMPLE_NONINDEX = 'ISPRESENT'


class Annotator(object):
    def __init__(self, configfile, ped, ex_cont_re):
        """
        :param ex_cont_re: List[regex]
         Exclude contigs matching any of these regexes from consideration
          Aka will _always_ return True for these regexes
        """
        self.config = set_config(configfile)
        self.readers = {}
        self.hgvs_parser = hgvs.parser.Parser()
        self.index = ped.get_index_samples()[0].individual_id
        self.exclude_contig_regex = ex_cont_re

    def __prep_vcf(self, path):
        if path not in self.readers:
            self.readers[path] = vcf.Reader(filename=path)
        return self.readers[path]

    def __prep_tabix(self, path):
        if path not in self.readers:
            self.readers[path] = Tabixfile(path)
        return self.readers[path]

    def __get_vcf_region(self, chrom, pos, reader, n_it=0):
        """
        Attempt to get vcf region
        :param chrom:
        :param pos:
        :param reader:
        :param n_it
        :return: iter or none
        """
        iter = None
        if n_it > 1:
            return iter

        try:
            # since pyvcf 0.6.8 regions are 0-indexed.
            if _is_vcf_version_at_least_0_6_8():
                iter = reader.fetch(chrom, pos-1, pos)
            else:
                iter = reader.fetch(chrom, pos, pos)
        except ValueError:
            # a little recursion
            if chrom.startswith('chr'):
                chr = chrom.split("chr")[1]
                iter = self.__get_vcf_region(chr, pos, reader, n_it=n_it+1)
            else:
                chr = "chr" + chrom
                iter = self.__get_vcf_region(chr, pos, reader, n_it=n_it+1)

        return iter

    def get_custom_annotations(self, record):
        annot = {}
        for x in self.config['vcf'].get('id', []):
            annot[x.get('colname', None)] = self.vcf_id(
                record,
                self.__prep_vcf(x.get('path', None))
            )

        for x in self.config['vcf'].get('count', []):
            annot[x.get('colname', None)] = self.vcf_count_occurences(
                record,
                self.__prep_vcf(x.get('path', None))
            )

        for x in self.config['vcf'].get('info', []):
            annot[x.get('colname', None)] = self.vcf_frequency(
                record,
                self.__prep_vcf(x.get('path', None)),
                x.get('ref_colname', None)
            )

        for x in self.config.get('tabix', []):
            annot.update(self.tabix_file(
                record,
                self.__prep_tabix(x.get('path', None)),
                x.get('colnames', []))
            )

        for x in self.config.get('bigwig', []):
            annot[x.get('colname', None)] = self.bigwig_directory(
                record, x.get('path', None)
            )

        return annot

    def vcf_count_occurences(self, record, reference_vcf):
        """
        Count the number of occurences of record in reference vcf
        :param record: vcf record
        :param reference_vcf: vcf reader of reference vcf
        :return: Count
        """
        ref_iter = self.__get_vcf_region(
            record.CHROM,
            record.POS,
            reference_vcf
        )

        if ref_iter is None:
            return 'unknown'
        count = 0
        db_vals = [(x.REF, [str(z) for z in x.ALT]) for x in ref_iter]
        for db_ref, db_alt in db_vals:
            if record.REF == db_ref and \
                            sorted([str(x) for x in record.ALT]) == sorted(db_alt):
                count += 1

        return collapse_values(count)

    def vcf_frequency(self, record, reference_vcf,
                      column_name, fallback='unknown'):
        """
        :param record:
        :param reference_vcf:
        :param column_name:
        :param fallback:
        :return:
        """
        ref_iter = self.__get_vcf_region(
            record.CHROM,
            record.POS,
            reference_vcf
        )

        if ref_iter is None:
            return fallback

        for ref_rec in ref_iter:
            col = reference_vcf.infos.get(column_name, None)
            if col is None:
                return fallback
            info_number = col.num
            val = self.__get_ref_record_column(
                record,
                ref_rec,
                column_name,
                info_number
            )
            if val is not None:
                return collapse_values(val)

        return fallback

    def __get_ref_record_column(self, rec, ref_rec, column,
                                info_number, fallback='unknown'):
        """
        Get the correct column value of a reference record
        :param rec: query record
        :param ref_rec: reference reocrd
        :param column: column name
        :param info_number: number of info field.
        Uses pyvcf-style numbers (i.e. A = -1, G = -2, R = -3)
        :return: value or None
        """

        s1, s2 = None, None

        if info_number is None:  # no info about alleles
            return collapse_values(ref_rec.INFO.get(column, fallback))
        elif info_number >= 0:  # no info about alleles
            return collapse_values(ref_rec.INFO.get(column, fallback))
        elif info_number == -1:  # one value per alt allele
            s1 = set([x.sequence for x in rec.ALT])
            s2 = set([x.sequence for x in ref_rec.ALT])
        elif info_number == -2:  # one value per possible genotype
            return NotImplementedError
        elif info_number == -3:  # one value per each possible allele
            s1 = set([x.sequence for x in rec.ALT]).union({rec.REF})
            s2 = set([x.sequence for x in ref_rec.ALT]).union({ref_rec.REF})

        inter = s1.intersection(s2)
        if len(inter) == 0:
            return fallback
        idxs = [list(s2).index(x) for x in inter]
        col = ref_rec.INFO.get(column, None)
        if col is None:
            return fallback
        if len(idxs) == 1:
            return collapse_values(col[idxs[0]])
        else:
            return collapse_values([col[x] for x in idxs])

    def vcf_id(self, record, reference_vcf, fallback='unknown'):
        iter = self.__get_vcf_region(record.CHROM, record.POS, reference_vcf)
        if iter is None:
            return fallback
        else:
            for rec in iter:
                if rec.ID is None:
                    return fallback
                else:
                    return rec.ID
            return fallback

    def tabix_file(self, record, reference_tabix,
                   columns=None, fallback='unknown'):
        """
        Retrieve columns from a tabix file
        (i.e. gzipped and tabixxed tsv file)
        :param record: query VCF record
        :param reference_tabix: instance of Tabixfile
        :param columns: names of columns to retrieve
        :param fallback: fallback value
        :return: dict of {column: value}
        """
        default = {x: fallback for x in columns}
        if record.CHROM.startswith('chr') and not any(
                [x.startswith('chr') for x in reference_tabix.contigs]
        ):
            chrom = record.CHROM.replace('chr', '')
        else:
            chrom = record.CHROM
        # adjust range (pysam's ranges are 0-based, half-open)
        start, end = record.POS - 1, record.POS

        # we ONLY consider the LAST line of any header field. This MUST be tab-delimited
        if not hasattr(self, "__tabix_header_{0}".format(reference_tabix)):
            setattr(self, "__tabix_header_{0}".format(reference_tabix),
                    [x for x in reference_tabix.header])
        tabix_header = getattr(
            self, "__tabix_header_{0}".format(reference_tabix)
        )[-1].strip().split("\t")

        try:
            ref_iter = reference_tabix.fetch(chrom, start=start, end=end)
        except ValueError:
            return default

        ref_piter = PeekableIterator(ref_iter)

        if ref_piter.peek() is None:
            return default

        for ref_rec in ref_piter:
            row = ref_rec.split('\t')
            assert len(row) == len(tabix_header), \
                "{0} {1}".format(len(row), len(tabix_header))
            # check POS, REF, ALT ~ only fetch value if all three are the same
            # with their record counterpart
            if str(record.POS) == row[1] and record.REF == row[2] and \
                    ','.join([str(x) for x in record.ALT]) == row[3]:
                for colname in columns:
                    default[colname] = collapse_values(row[tabix_header.index(colname)])
            else:
                if ref_piter.peek() is None:
                    return default

        return default

    def bigwig_directory(self, record, reference_bw_dir, fallback='unknown'):
        """
        Query value from a folder containing bigwig tracks

        bigwig names must be *IDENTICAL* to the contig names in your record
        :param record: query record
        :param reference_bw_dir: folder containing bigwig files
        :param fallback: fallback value
        :return: value
        """
        bw_name = os.path.join(reference_bw_dir, record.CHROM + '.bw')
        # query bw files using pyBigWig ~ TAKE COORDINATE SETTING INTO ACCOUNT
        bw = pyBigWig.open(bw_name)
        qres = bw.intervals(record.CHROM, record.POS - 1, record.POS)
        if qres:
            return collapse_values(qres[0][-1])

        return fallback

    def _get_record_values(self, rec, fallback='unknown'):

        d = {}
        for column in COLS_RECORD:
            try:
                d[column] = collapse_values(rec.INFO[column][0])
            except:
                d[column] = fallback

        return d

    def _get_sample_columns(self, rec, fallback='unknown'):
        returnval = {}
        for sample in rec.samples:
            # calculate AD first
            ad = getattr(sample.data, 'AD', None)
            altperc = 0.0
            dpref = 0
            dpalt = 0
            # only if Allele Depth is available and not both to zero (sum>1)
            if ad and sum(ad):
                altperc = round(ad[1] * 100 / float(sum(ad)), 3)
                dpref = ad[0]
                dpalt = ad[1]
            else:
                ad = 'unknown'
            ad_data = {'AD': ad, 'ALTPERC': altperc, 'DP': dpalt + dpref,
                       'DPALT': dpalt, 'DPREF': dpref}

            for ps_colname in COLS_PER_SAMPLE:
                if ps_colname in ('GT', 'GQ'):
                    col_value = getattr(sample.data, ps_colname, fallback)
                    if col_value is None:
                        value = fallback
                    else:
                        value = col_value
                elif ps_colname in ('AD', 'ALTPERC', 'DPALT', 'DPREF', 'DP'):
                    value = ad_data[ps_colname]
                elif ps_colname == 'PLALT':
                    pl = getattr(sample.data, 'PL', fallback)
                    if pl is not None:
                        value = pl[1]
                    else:
                        value = fallback
                else:
                    value = fallback
                returnval[ps_colname + '_' + sample.sample] = value

        return returnval

    def _get_sample_nonindex_columns(self, rec, line):
        # then per non-index sample columns
        returnval = {}
        for sample in rec.samples:
            sname = sample.sample
            if sname != self.index:
                altperc = line['ALTPERC_' + sname]
                plalt = line['PLALT_' + sname]
                colname = COLS_PER_SAMPLE_NONINDEX + '_' + sname

                # default value is 1
                value = '1'
                gt = getattr(sample.data, 'GT')
                if gt is not None and (not gt.startswith('0') or not
                   gt.endswith('0')):
                    value = '6'
                elif altperc > 10:
                    value = '5'
                elif 0 < altperc <= 10:
                    value = '4'
                elif plalt < 30 or plalt == 'unknown':
                    value = '3'
                elif plalt < 60:
                    value = '2'

                returnval[colname] = value
        return returnval

    def is_refseq(self, transcript):
        """
        Checks whether transcript is a refseq transcripts
        """
        contig = transcript.CHROM
        if any([r.match(contig) for r in self.exclude_contig_regex]):
            return True
        feature = transcript.INFO['Feature']
        if feature.startswith('ENS'):
            return False
        elif feature.startswith('CCDS'):
            return False
        elif feature.startswith('XM'):
            return True
        elif feature.startswith('NR'):
            return False
        elif feature.startswith('NM'):
            return True
        else:
            return False

    def _distance_to_splice(self, transcript):
        """
        Get distance to splice site of an intron transcript
        Note: this assumes we have an intronic variant!
        Do note attempt this with other variant types!
        :param transcript: vcf record with Feature info tag
        :return: distance to splice site (integer) or None if unknown
        """

        description = transcript.INFO['HGVSc']
        if description is None or description == '':
            return None

        try:
            variant = self.hgvs_parser.parse_hgvs_variant(description)
        except HGVSParseError:
            sys.stderr.write("WARNING! HGVS notation {0} "
                             "seems to be invalid\n".format(description))
            return None

        base_offset = variant.posedit.pos.start

        base = base_offset.base
        offset = base_offset.offset

        # intronic
        if base > 0:
            return offset
        # upstream of start codon
        elif base < 0 and offset == 0:
            return base
        # downstream of stop codon
        elif base == 0 and offset > 0:
            return offset
        else:
            return None

    def _get_distance_to_splice(self, transcript, colname):
        """
        Wrapper function for formatting dsp into common dictionary object
        :param transcript: vcf record
        :param colname: column name to be used
        :return: dict of {"colname": dsp_value}
        """

        # select for intronic
        if "intron_variant" in transcript.INFO["Consequence"]:
            dsp = self._distance_to_splice(transcript)
            if dsp is not None:
                return {colname: dsp}
            else:
                return {colname: 'unknown'}

        else:
            return {colname: 'unknown'}
