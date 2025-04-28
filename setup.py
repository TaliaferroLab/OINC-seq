#from distutils.core import setup
#add tests?
from setuptools import setup, find_packages
setup(name = 'pigpen',
description = 'Pipeline for the Identification of Guanosine Positions Erroneously Notated',
author = 'Matthew Taliaferro',
author_email = 'taliaferrojm@gmail.com',
url = 'https://github.com/TaliaferroLab/OINC-seq',
version = '0.0.7',
packages = find_packages(where = './src', exclude = ['workflow', 'testdata']),
package_dir = {'':'src'},
entry_points = {'console_scripts': ['pigpen = pigpen.runpigpen:main', 'bacon = pigpen.bacon_glm:main', 'alignAndQuant = pigpen.alignAndQuant:main']})
