from setuptools import setup, find_packages

setup(
    name="smorfinder",
    version='1.0.0',
    description='A command line tool to identify and annotate small proteins in genomes and metagenomes.',
    url='https://github.com/bhattlab/SmORFinder',
    author="Matt Durrant",
    author_email="mdurrant@stanford.edu",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'click==7.0',
        'biopython==1.76',
        'tensorflow==2.3.1',
        'wget==3.2'
    ],
    zip_safe=False,
    entry_points = {
        'console_scripts': [
            'smorf = smorfinder.main:cli'
        ],
}
)
