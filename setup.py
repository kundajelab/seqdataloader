from setuptools import setup,find_packages

config = {
    'include_package_data': True,
    'author': 'Anna Shcherbina',
    'author_email': 'annashch@stanford.edu',
    'url': 'https://github.com/kundajelab/seqdataloader',
    'description': 'Generate genome-wide classification and regression labels for DNA accessibility data.',
    'version': '0.113',
    'packages': ['seqdataloader'],
    'setup_requires': [],
    'install_requires': ['numpy','pandas','cython','deeptools','pybedtools','pyBigWig'],
    'scripts': [],
    'entry_points': {'console_scripts': ['genomewide_labels=seqdataloader.__init__:main']},
    'name': 'seqdataloader'
}

if __name__== '__main__':
    setup(**config)
