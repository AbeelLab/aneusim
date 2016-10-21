from setuptools import setup

setup(
    name='aneuploidgen',
    version='0.1',
    py_modules=['aneuploidgen'],

    # Metadata
    description='Helper script to generate synthetic aneuploid genomes.',
    author='Lucas van Dijk',
    author_email='info@lucasvandijk.nl',
    license='MIT',
    url='https://bitbucket.org/tudelft-bioinformatics/aneugen',

    # Dependencies
    install_requires=[
        'scikit-bio>=0.5.0'
    ],
    setup_requires=['pytest-runner'],
    tests_require=[
        'pytest'
    ],

    # Entry points
    entry_points={
        'console_scripts': [
            'aneuploidgen = aneuploidgen:main'
        ]
    }
)
