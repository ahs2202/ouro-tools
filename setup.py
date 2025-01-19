"""Setup script for Ouro-Tools """

import os.path
from setuptools import setup, find_packages

# The directory containing this file
HERE = os.path.abspath(os.path.dirname(__file__))

# The text of the README file
with open(os.path.join(HERE, "README.md")) as fid:
    README = fid.read()

setup(
    name="ourotools",
    version="0.2.5",
    author="Hyunsu An",
    author_email="ahs2202@gm.gist.ac.kr",
    description="A comprehensive toolkit for quality control and analysis of single-cell long-read RNA-seq data",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/ahs2202/ouro-tools",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    python_requires=">=3.11, <4",
    install_requires=[
        "bitarray>=2.5.1",
        "scipy>=1.9.1",
        "tqdm>=4.64.1",
        "nest-asyncio>=1.5.6",
        "joblib>=1.2.0",
        "pandas>=1.5.2",
        "intervaltree>=3.1.0",
        "matplotlib>=3.5.2",
        "plotly>=5.18.0",
        "regex>=2.5.135",
        "owlready2>=0.46",
        "scanpy>=1.10.2",
        "h5py>=3.8.0",
        "mappy>=2.24",
        "pysam>=0.18.0",
        "pyBigWig>=0.3.22",
    ],
    entry_points={
        "console_scripts": [
            "ourotools=ourotools.core.core:ourotools",
        ]
    },
)
