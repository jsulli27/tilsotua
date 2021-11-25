from __future__ import absolute_import, division, print_function

# Standard imports
import sys
import glob, os
from setuptools import setup, find_packages


# Begin setup
setup_keywords = dict()
setup_keywords['name'] = 'tilsotua'
setup_keywords['description'] = 'Inverse AUTOSLIT for LRIS masks'
setup_keywords['author'] = 'Jessica Sullivan et al'
setup_keywords['author_email'] = 'jsulli27@nd.edu'
setup_keywords['license'] = ''
setup_keywords['url'] = 'https://github.com/jsulli27/tilsotua'
setup_keywords['version'] = '0.0.dev0'
# Use README.rst as long_description.
setup_keywords['long_description'] = ''
if os.path.exists('README.md'):
    with open('README.md') as readme:
        setup_keywords['long_description'] = readme.read()
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['requires'] = ['Python (>3.7.0)']
setup_keywords['install_requires'] = ['astropy>=4.3.post1', 'matplotlib>=3.4.2','shapely',
    'panstamps>=0.6.1', 'pytest-shutil>=1.7.0', 'numpy>=1.21.1','astroquery>=0.4.3','regions>=0.5','pyyaml>=5.3.1']
setup_keywords['zip_safe'] = False
setup_keywords['use_2to3'] = False
setup_keywords['packages'] = find_packages()
setup_keywords['setup_requires'] = ['pytest-runner']
setup_keywords['tests_require'] = ['pytest']

if os.path.isdir('bin'):
    setup_keywords['scripts'] = [fname for fname in glob.glob(os.path.join('bin', '*'))
                                 if not os.path.basename(fname).endswith('.rst')]

setup(**setup_keywords)