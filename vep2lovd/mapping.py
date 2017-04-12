"""
vep2lovd.mapping
~~~~~~~~~~~~~
:copyright: (c) 2016 Sander Bollen
:copyright: (c) 2016 Leiden University Medical Center
:license: MIT
"""
__author__ = 'Sander Bollen'


class GVSFunction(object):
    """
    Class for mapping VEP consequences to SSEQ consequences
    """
    def __init__(self):
        self._map = {'3_prime_UTR_variant ': 'utr-3',
                     '5_prime_UTR_variant': 'utr-5',
                     'coding_sequence_variant': 'coding',
                     'coding_sequence_variant&3_prime_UTR_variant': 'codingComplex',
                     'coding_sequence_variant&5_prime_UTR_variant': 'codingComplex',
                     'frameshift_variant&stop_lost': 'frameshift',
                     'frameshift_variant': 'frameshift',
                     'frameshift_variant&splice_region_variant': "frameshift",
                     'incomplete_terminal_codon_variant&coding_sequence_variant': 'coding-near-splice',
                     'inframe_deletion': 'coding',
                     'inframe_deletion&splice_region_variant': 'coding-near-splice',
                     'inframe_deletion&stop_retained_variant': 'coding',
                     'inframe_insertion': 'coding',
                     'inframe_insertion&splice_region_variant': 'coding-near-splice',
                     'initiator_codon_variant': 'coding',
                     'initiator_codon_variant&splice_region_variant': 'coding-near-splice',
                     'intron_variant': 'intron',
                     'missense_variant': 'missense',
                     'missense_variant&splice_region_variant': 'missense-near-splice',
                     'splice_acceptor_variant': 'splice-3',
                     'splice_acceptor_variant&5_prime_UTR_variant ': 'utr-5',
                     'splice_acceptor_variant&coding_sequence_variant': 'splice-3',
                     'splice_acceptor_variant&coding_sequence_variant&intron_variant': 'splice-3',
                     'splice_acceptor_variant&intron_variant': 'splice-3',
                     'splice_acceptor_variant&splice_donor_variant&coding_sequence_variant&intron_variant':
                         'codingComplex-near-splice',
                     'splice_acceptor_variant&splice_donor_variant&inframe_deletion&intron_variant':
                         'codingComplex-near-splice',
                     'splice_donor_variant': 'splice-5',
                     'splice_donor_variant&5_prime_UTR_variant&intron_variant': 'utr-5',
                     'splice_donor_variant&coding_sequence_variant': 'splice-5',
                     'splice_donor_variant&coding_sequence_variant&intron_variant': 'splice-5',
                     'splice_donor_variant&intron_variant': 'splice-5',
                     'splice_region_variant&3_prime_UTR_variant': 'utr-3',
                     'splice_region_variant&5_prime_UTR_variant': 'utr-5',
                     'splice_region_variant&intron_variant': 'splice',
                     'splice_region_variant&synonymous_variant': 'coding-synonymous-near-splice',
                     'stop_gained': 'stop-gained',
                     'stop_gained&inframe_deletion': 'stop-gained',
                     'stop_gained&inframe_insertion': 'stop-gained',
                     'stop_gained&inframe_insertion&splice_region_variant': 'stop-gained',
                     'stop_gained&splice_region_variant': 'stop-gained',
                     'stop_lost': 'stop-lost',
                     'stop_lost&inframe_deletion': 'stop-lost',
                     'stop_retained_variant': 'coding-synonymous',
                     'synonymous_variant': 'coding-synonymous',
                     'upstream_gene_variant': 'utr-5',
                     'downstream_gene_variant': "utr-3",
                     '3_prime_UTR_variant': 'utr-3',
                     'splice_donor_variant&3_prime_UTR_variant&intron_variant': 'utr-3'}

    def map(self, consequence):
        try:
            return self._map[consequence]
        except KeyError:
            return "unknown"


DEFAULT_GATK_MAPPING = {
    "UG": "UG",
    "HC": "HC",
    "UG,HC": "UG,HC"
}
