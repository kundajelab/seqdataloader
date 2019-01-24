from setuptools import setup,find_packages

config = {
    'include_package_data': True,
    'description': 'Generate genome-wide classification and regression labels for DNA accessibility data.',
    'version': '0.21',
    'setup_requires': [],
    'install_requires': ['numpy','pandas','cython','deeptools','pybedtools','pyBigWig'],
    'scripts': [],
    'entry_points': {'console_scripts': ['genomewide_labels = genomewide_labels.py:main']},
    'name': 'seqdataloader'
}

if __name__== '__main__':
    setup(**config)
