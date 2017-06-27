"""
vep2lovd.filter
~~~~~~~~~~~~~
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: MIT
"""
__author__ = 'Sander Bollen'


class MultiFilter(object):
    """
    Class for filtering out multiple transcripts at once
    """

    def __init__(self, annotator):
        self.annotator = annotator

    def filter_xm(self, transcripts):
        refseq_transcripts = [x for x in transcripts if
                              self.annotator.is_refseq(x)]
        nm_transcript = []
        for tr in refseq_transcripts:
            if tr.INFO['Feature'].startswith("NM"):
                nm_transcript.append(tr)
        if len(nm_transcript) > 0:
            return nm_transcript
        else:
            return refseq_transcripts

    def filter_same_transcripts(self, transcripts):
        """
        This function filters out transcripts which are identical
        """
        features = []
        trs = []
        for tr in transcripts:
            if tr.INFO['Feature'] not in features:
                features.append(tr.INFO['Feature'])
                trs.append(tr)
        return trs


class SingleFilter(object):
    """
    Class for filtering based on a single transcript
    """
    def __init__(self, annotator, filter_consequences=None,
                 filter_contains_consequences=None):
        self.annotator = annotator
        self.filter_cons = filter_consequences
        self.filter_contains_cons = filter_contains_consequences

    def filter_consequences(self, transcript):
        """
        This function returns false if transcript is in FILTER_CONSEQUENCES
        ignore for excluded contigs
        """
        contig = transcript.CHROM
        if any([r.match(contig) for r in self.annotator.exclude_contig_regex]):
            return True
        if transcript.INFO["Consequence"] in self.filter_cons:
            return False
        else:
            return True

    def filter_contains_consequences(self, transcript):
        """
        This function retiurns false if transcript has consequence
        matching anything in FILTER_CONTAINS_CONSEQUENCES
        ignore for exclude contigs
        """
        val = True
        contig = transcript.CHROM
        if any([r.match(contig) for r in self.annotator.exclude_contig_regex]):
            return True
        for filter in self.filter_contains_cons:
            if filter in transcript.INFO["Consequence"]:
                val = False

        return val

    def filter_distance_to_splice(self, transcript, maximum):
        """
        Check whether a transcript has a distance to splice greater than max
        Non-intronic transcripts will always return False
        :param transcript: vcf record with assumed Consequence and
        Feature info fields
        :param maximum: integer, max allowed dsp
        :return: Boolean
        """

        # check intronic
        if "intron_variant" in transcript.INFO['Consequence']:

            try:
                dsp = abs(self.annotator._distance_to_splice(transcript))
            except TypeError:
                return False

            if dsp is not None:
                return dsp > maximum
            else:
                return False
        else:
            return False
