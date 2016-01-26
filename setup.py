#!/bin/env python

from distutils.core import setup

setup(name='hts_scripts',
      version='0.5.0',
      description='Scripts for processing NGS data',
      author='Rob Carter',
      author_email='robert.carter@stjude.org',
      #packages=['hlatyp'],
      #package_data={'hlatyper': ['data/*']},
      scripts=['scripts/sample_bam.py']
     )
