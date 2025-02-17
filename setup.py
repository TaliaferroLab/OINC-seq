#from distutils.core import setup
from setuptools import setup, find_packages
setup(name = 'pigpen',
description = 'Pipeline for the Identification of Guanosine Positions Erroneously Notated',
author = 'Matthew Taliaferro',
author_email = 'taliaferrojm@gmail.com',
url = 'https://github.com/TaliaferroLab/OINC-seq',
version = '0.0.2',
packages = find_packages('.', exclude = ['workflow', 'testdata']),
scripts = ['ExtractUMI.py', 'alignAndQuant.py', 'alignUMIquant.py', 'assignreads.py', 'assignreads_salmon.py', 'assignreads_salmon_ensembl.py', 'bacon_glm.py', 'conversionsPerGene.py', 'filterbam.py', 'getmismatches.py', 'getmismatches_MPRA.py', 'maskpositions.py', 'parsebamreadcount.py', 'pigpen.py', 'snps.py'])
