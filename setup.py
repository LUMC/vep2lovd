__author__ = 'ahbbollen'

from setuptools import setup

setup(name="vep2lovd",
      version="0.1",
      description="Vep2lovd parser",
      author="Sander Bollen",
      author_email="a.h.b.bollen@lumc.nl",
      packages=['vep2lovd'],
      zip_safe=False,
      install_requires=['numpy', 'pyvcf', 'hgvs==0.3.7', 'pysam', 'bx-python==0.7.1', "six"])