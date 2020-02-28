#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

# TODO - Should be install via conda or here?? biopython?
# TODO - pyyaml it's the conda install pyyaml and not "yaml"!!!
# TODO - edlib & send2trash is fucked up...no version
requirements = ['Click>=6.0', 'matplotlib>=3.1.0', 'seaborn>=0.9.0','pyyaml>=5.1.2', 'plotly>=4.3.0'
                'numpy>=1.12.1', 'pandas>=0.24.2', 'edlib', 'biopython==1.74', 'scipy>=1.2.1', 'statsmodels',
                'jinja2']

setup_requirements = ['pytest-runner']

test_requirements = ['pytest']

setup(
    author="Ido Amit",
    author_email='idoamit98@gmail.com',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="CRISPECTOR - Genome Editing Analysis Tool",
    entry_points={
        'console_scripts': [
            'crispector=crispector.arg_parser:main',
        ],
    },
    install_requires=requirements,
    license="GNU General Public License v3",
    include_package_data=True,
    keywords='crispector',
    name='crispector',
    packages=find_packages(include=['crispector']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/iamit87/crispector',
    version='0.1.0',
    zip_safe=False,
)
