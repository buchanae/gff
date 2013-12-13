from distutils.core import setup
import os

_this_dir = os.path.dirname(__file__)
README_path = os.path.join(_this_dir, 'README.md')
README = open(README_path).read()


setup(
    name='gff',
    description='Tools for reading and working with GFF v3 records',
    long_description=README,
    version='2.0.0',
    author='Alex Buchanan',
    author_email='buchanae@gmail.com',
    license='MIT',
    py_modules=[
        'gff.formatter',
        'gff.parser',
        'gff.record',
        'gff.tree',
    ]
)
