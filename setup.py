#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup

requirements = ['Click>=7.0', 'matplotlib >= 3.1.2', 'seaborn>=0.9.0','pyyaml>=5.1.2', 'plotly>=4.3.0'
                'numpy>=1.12.1', 'pandas>=0.24.2', 'biopython>=1.74', 'scipy>=1.2.1', 'statsmodels',
                'jinja2']

setup(
    author="Ido Amit",
    author_email='idoamit98@gmail.com',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: Other/Proprietary License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.7',

    ],
    description="CRISPECTOR - Genome Editing Analysis Tool",
    entry_points={
        'console_scripts': [
            'crispector = crispector.cli:main',
        ],
    },
    python_requires='>=3.7',
    install_requires=requirements,
    include_package_data=True,
    keywords='crispector',
    name='crispector',
    package_dir={'crispector': 'crispector'},
    packages=['crispector', 'crispector.algorithm', 'crispector.config', 'crispector.input_processing',
              'crispector.modifications', 'crispector.report', 'crispector.utils'],
    url='https://github.com/YakhiniGroup/crispector',
    version='1.0.2b7',
    zip_safe=False,
)