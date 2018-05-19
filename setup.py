# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 10:39:53 2017

@author: Daniel Uher

Setup file for refsig package
"""

from setuptools import setup

setup(name='refsig',
      version='1.0.0',
      description='Package for obtaining the referential signal from a set of iEEG data',
      author='DandaPanda',
      author_email='daniel-uher@hotmail.com',
      url='https://github.com/DandaPanda/refsig',
      license='MIT',
      packages=['refsig'],
      keywords='reference iEEG EEG referential signal ICA',
      zip_safe=False)