from setuptools import setup, find_packages

setup(
    name="deepsmorfnet",
    version='0.0.12_dev',
    description='A command line tool to identify and annotate small proteins in metagenomic sequencing datasets.',
    url='https://github.com/bhattlab/DeepSmORFNET',
    author="Matt Durrant",
    author_email="mdurrant@stanford.edu",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'click==7.0',
        'biopython==1.76',
        'tensorflow==2.1.0',
        'wget==3.2'
    ],
    zip_safe=False,
    entry_points = {
        'console_scripts': [
            'dsn = deepsmorfnet.main:cli'
        ],
}
)
