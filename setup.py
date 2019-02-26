from setuptools import setup,find_packages
config = {
    'include_package_data': True,
    'author': 'Anna Shcherbina',
    'author_email': 'annashch@stanford.edu',
    'url': 'https://github.com/kundajelab/seqdataloader',
    'description': 'Generate genome-wide classification and regression labels for DNA accessibility data.',
    'version': '0.120',
    'packages': ['seqdataloader'],
    'setup_requires': [],
    'install_requires': ['numpy>=1.15','pandas>=0.23.4','cython>=0.27.3','deeptools>=3.0.1','pybedtools>=0.7','pyBigWig>=0.3.2'],
    'scripts': [],
    'entry_points': {'console_scripts': ['genomewide_labels=seqdataloader.__init__:main']},
    'name': 'seqdataloader'
}

if __name__== '__main__':
    setup(**config)
