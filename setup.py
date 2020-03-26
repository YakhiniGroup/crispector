#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

requirements = ['Click>=6.0', 'matplotlib>=3.1.0', 'seaborn>=0.9.0','pyyaml>=5.1.2', 'plotly>=4.3.0'
                'numpy>=1.12.1', 'pandas>=0.24.2', 'biopython==1.74', 'scipy>=1.2.1', 'statsmodels',
                'jinja2', 'binascii', 'gzip']

setup(
    author="Ido Amit",
    author_email='idoamit98@gmail.com',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: Other/Proprietary License',
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
            'crispector = crispector.arg_parser:cli',
        ],
    },
    install_requires=requirements,
    include_package_data=True,
    keywords='crispector',
    name='crispector',
    packages=find_packages(include=['crispector']),
    url='https://github.com/YakhiniGroup/crispector',
    version='1.0.1b',
    zip_safe=False,
)
