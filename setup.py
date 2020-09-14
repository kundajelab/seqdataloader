from setuptools import setup,find_packages
config = {
    'include_package_data': True,
    'author': 'Anna Shcherbina',
    'author_email': 'annashch@stanford.edu',
    'url': 'https://github.com/kundajelab/seqdataloader',
    'description': 'Generate genome-wide classification and regression labels for DNA accessibility data.',
    'version': '1.0',
    'packages': ['seqdataloader'],
    'setup_requires': [],
    'install_requires': ['numpy>=1.15','pandas>=0.23.4','cython>=0.27.3','deeptools>=3.0.1','pybedtools>=0.7','pyBigWig>=0.3.7', 'pyfaidx','tiledb>=0.4.4'],
    'scripts': [],
    'entry_points': {'console_scripts': ['genomewide_labels=seqdataloader.labelgen.__init__:main',
                                         'db_ingest=seqdataloader.dbingest.__init__:main',
                                         'db_ingest_single_threaded=seqdataloader.dbingest_single_threaded.__init__:main',
                                         'seqdataloader_get_outliers=seqdataloader.helpers.get_outliers:main']},
    'name': 'seqdataloader'
}

if __name__== '__main__':
    setup(**config)
