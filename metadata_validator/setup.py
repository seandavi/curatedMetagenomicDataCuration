"""
Setup script for the Metadata Validator Framework

DEPRECATED: This file is kept for backward compatibility only.
Please use pyproject.toml with uv or pip instead:
    uv pip install -e .
    pip install -e .

See UV_GUIDE.md for more information about using uv (much faster!).
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="metadata-validator",
    version="0.1.0",
    author="Curated Metagenomic Data Curation Team",
    description="A comprehensive Python framework for validating, harmonizing, and enriching sample metadata",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "pandas>=1.3.0",
        "numpy>=1.21.0",
    ],
    extras_require={
        "ontology": [
            "pronto>=2.5.0",
            "requests>=2.26.0",
        ],
        "fuzzy": [
            "python-Levenshtein>=0.12.2",
            "fuzzywuzzy>=0.18.0",
        ],
        "dev": [
            "pytest>=6.0",
            "black>=21.0",
            "flake8>=3.9",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
)
