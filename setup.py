from setuptools import setup, find_packages

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

from cgevirus.version import __version__

setup(
    name='cgevirus',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version=__version__,
    packages=find_packages(),
    data_files=[],
    include_package_data=True,
    url='https://https://github.com/MBHallgren/cgevirus',
    license='',
    install_requires=(),
    author='Malte B. Hallgren',
    scripts=['bin/cgevirus'],
    author_email='malhal@food.dtu.dk',
    description='cgevirus - K-mer Gene Typer',
)