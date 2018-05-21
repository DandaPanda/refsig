# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 10:39:53 2017

@author: Daniel Uher

Setup file for refsig package
"""

from setuptools import setup
from codecs import open
from os import path

def readme():
    with open('README.rst') as f:
        return f.read()

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()
    
    
setup(name='refsig',
      version='1.0.4',
      description='Package for obtaining the referential signal from a set of iEEG data',
      long_description = long_description,
      author='Daniel Uher in affiliation with FNUSA-ICRC Brno and BUT Brno',
      author_email='daniel-uher@hotmail.com',
      url='https://github.com/DandaPanda/refsig',
      license='MIT',
      packages=['refsig'],
      keywords='reference iEEG EEG referential signal ICA',
      zip_safe=False)