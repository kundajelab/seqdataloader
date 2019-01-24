from setuptools import setup,find_packages

config = {
    'include_package_data': True,
    'description': 'Generate genome-wide classification and regression labels for DNA accessibility data.',
    'version': '0.1',
    'packages': ['seqdataloader'],
    'setup_requires': [],
    'install_requires': ['numpy','pybedtools','pyBigWig','pandas','multiprocessing'],
    'scripts': [],
    'entry_points': {'console_scripts': ['genomewide_labels = seqdataloader.genomewide_labels.py:main']},
    'name': 'seqdataloader'
}

if __name__== '__main__':
    setup(**config)
