VEP2LOVD
=========

This is a package for annotating, filtering, and converting VEP-annotated VCFs to LOVD format.

LOVD is transcript-based rather than variant-based, so we will have to make some implicit conversions.
We therefore expect the input to this module to be pre-filtered; i.e. removing all records for which no VEP annotation exists (e.g. Mitochondrial DNA)

Installation
------------

To install vep2lovd, we recommend using a virtual environment.

`numpy` _must_ be installed before running `pip install -r requirements.txt`.
After installing `numpy`, install all other dependencies with `pip install -r requirements.txt`
Finally, install vep2lovd with `python setup.py develop`


How to use
----------

The module is split in several submodules:

1. `annotate` : contains an `Annotator` object to annotate records with additional info
2. `filter`: contains a `SingleFilter` object for filtering single transcripts, and a `MultiFilter` for group-based filtering
3. `mapping`: contains a `GVSFunction` object to map VEP consequences to GVS function names
4. `utils`: contains several helper functions and classes
5. `parser`: module for converting to LOVD format, including `LOVDFile` object. Generally, this will likely be your entry-point.


Convert VEP-annotated VCF to LOVD
----------------------------------

The following snippet of code will convert a VEP-annotated VCF to LOVD

```python
from vep2lovd import LOVDFile

ped_file = "/path/to/ped"
lovd_settings = "/path/to/settings.json"
vcf_file = "/path/to/input.vcf" # can be gzipped

lovd = LOVDFile(vcf_file, ped_file, lovd_settings)

with open("/path/to/output.data.lovd", "wb") as ohandle:
    for line in lovd:
        ohandle.write(line + "\n")

```

The constructor requires at least those three variables.

* `ped_file` is a pedigree file in PED format. See http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
* `settings_json` is a json containing paths to all the files used for annotation . 


The following parameters to the constructor are optional, but change the converter's behaviour

* `filter_consequences`: A list of consequences to filter out. This is a *soft* filter; only when the consequence length for the given transcript is 1, and it matches any of this list, the transcript will be filtered away.
* `filter_contains_consequences`: A list of consequences to filter out. This is a *hard* filter; when any consequence of the transcript matches anything in this list, the transcript will be removed
* `filter_dsp`: Integer with max distance to splice site. Default is 50. Set to 0 to disable.
* `filter_xm`: Boolean whether to filter out XM transcripts if an NM transcript for that variant also exists. Default is True
* `filter_refseq`: Boolean whether to filter out non-refseq transcript. Default is True
* `filter_identical`: Boolean whether to filter out identical transcripts. Default is True
