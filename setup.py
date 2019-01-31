from setuptools import setup,find_packages

config = {
    'include_package_data': True,
    'author': 'Anna Shcherbina',
    'author_email': 'annashch@stanford.edu',
    'url': 'https://github.com/kundajelab/seqdataloader',
    'description': 'Generate genome-wide classification and regression labels for DNA accessibility data.',
    'version': '0.11',
    'packages': ['genomewide_labels'],
    'setup_requires': [],
    'install_requires': ['numpy','pandas','cython','deeptools','pybedtools','pyBigWig'],
    'scripts': [],
    'entry_points': {'console_scripts': ['genomewide_labels=genomewide_labels.__init__:main']},
    'name': 'seqdataloader'
}

if __name__== '__main__':
    setup(**config)
