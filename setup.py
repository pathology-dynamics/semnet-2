
'''
setup.py file necessary to install SemNet package. To install, follow
instructions in readme.
'''

from setuptools import setup

setup(name='semnet',
	version='2.0.0',
	description='A package for working with Semantic Medline data',
	url='https://github.com/pathology-dynamics/semnet-2',
	author='Anna Kirkpatrick',
	author_email='akirkpatrick3@gatech.edu',
	packages=['semnet'],
	install_requires=['numpy>=1.20.3']
)
